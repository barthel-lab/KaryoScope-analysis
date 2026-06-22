"""Resolve overlapping per-featureset annotations to a single feature.

The ``overlay-annotations`` tool collapses several featureset annotation tracks
(region, repeat, subtelomeric, ...) into one feature per position. *How* they
collapse is defined by an explicit **resolution spec** — a file, not buried
``if/then`` (decision M2). To avoid the combinatorial blow-up of enumerating every
feature tuple, a spec is compact:

* ``precedence`` — the featuresets, in **default-winner** order. The first
  featureset wins wherever no rule matches.
* ``rules`` — an ordered list of exceptions. The **first** rule whose ``when``
  pattern matches a segment decides its output; ``when`` matches feature **names**,
  **lists**, or **hierarchy classes** (``@arm`` = ``arm`` and its descendants).

``emit`` forms:

* ``<featureset>`` — pass through that featureset's feature at the segment.
* ``{literal: X}`` — force feature ``X``.
* ``composite`` — join the matched ``when`` featuresets, in precedence order
  (default separator ``_``), e.g. ``DJ_TAR1``.
* ``{composite: all, sep: ":"}`` — join *all* precedence featuresets (the basic
  overlay/“composite” mode, e.g. ``region:repeat``).
* ``{composite: [chromosome, subtelomeric], sep: ":"}`` — join an *explicit* list of
  featuresets, in the given order (e.g. ``chr13:canonical_telomere``); lets a rule tag a
  resolved structural feature with another layer such as chromosome.

Specs are validated structurally (jsonschema) and semantically against the database
:class:`~karyoscope_analysis.core.feature_vocab.FeatureHierarchy` (every featureset,
feature name, and ``@class`` must exist), so an ill-formed spec fails before a run.
``when`` keys must be written in ``precedence`` order (a readability lint, M2).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

import jsonschema

from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

#: Reserved ``emit`` word; no featureset may be named this.
COMPOSITE = "composite"
_DEFAULT_SEP = {"when": "_", "all": ":", "explicit": ":"}

#: jsonschema for the *structure* of a spec (semantics are checked in the loader).
SPEC_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["name", "precedence"],
    "additionalProperties": False,
    "properties": {
        "name": {"type": "string"},
        "precedence": {"type": "array", "items": {"type": "string"}, "minItems": 1},
        "rules": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["emit"],
                "additionalProperties": False,
                "properties": {
                    "when": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {
                            "oneOf": [
                                {"type": "string"},
                                {"type": "array", "items": {"type": "string"}, "minItems": 1},
                            ]
                        },
                    },
                    "emit": {
                        "oneOf": [
                            {"type": "string"},
                            {
                                "type": "object",
                                "additionalProperties": False,
                                "properties": {
                                    "literal": {"type": "string"},
                                    "composite": {
                                        "oneOf": [
                                            {"enum": ["when", "all"]},
                                            {
                                                "type": "array",
                                                "items": {"type": "string"},
                                                "minItems": 1,
                                            },
                                        ]
                                    },
                                    "sep": {"type": "string"},
                                },
                            },
                        ]
                    },
                },
            },
        },
    },
}


class SpecError(ValueError):
    """A resolution spec is structurally or semantically invalid."""


@dataclass(frozen=True)
class Emit:
    """How a matched rule produces its output feature."""

    kind: Literal["featureset", "literal", "composite"]
    value: str | None = None  # featureset name (featureset) or literal value (literal)
    of: Literal["when", "all", "explicit"] = "when"  # composite scope
    sep: str = "_"
    featuresets: tuple[str, ...] = ()  # explicit composite: the featuresets to join, in order


@dataclass(frozen=True)
class Rule:
    """An ordered exception: ``when`` (per-featureset allowed features) → ``emit``."""

    when: Mapping[str, frozenset[str]]
    emit: Emit

    def matches(self, segment: Mapping[str, str]) -> bool:
        return all(segment.get(fs) in allowed for fs, allowed in self.when.items())


@dataclass(frozen=True)
class ResolutionSpec:
    """A loaded, validated overlay resolution spec."""

    name: str
    precedence: tuple[str, ...]
    rules: tuple[Rule, ...]

    def resolve(self, segment: Mapping[str, str]) -> str:
        """Resolve one refined segment (``{featureset: feature}``) to a single feature."""
        for rule in self.rules:
            if rule.matches(segment):
                return self._apply(rule, segment)
        return segment[self.precedence[0]]  # precedence default winner

    def _apply(self, rule: Rule, segment: Mapping[str, str]) -> str:
        emit = rule.emit
        if emit.kind == "featureset":
            assert emit.value is not None
            return segment[emit.value]
        if emit.kind == "literal":
            assert emit.value is not None
            return emit.value
        if emit.of == "explicit":
            featuresets = emit.featuresets
        elif emit.of == "all":
            featuresets = self.precedence
        else:  # "when"
            featuresets = tuple(fs for fs in self.precedence if fs in rule.when)
        return emit.sep.join(segment[fs] for fs in featuresets)


# --------------------------------------------------------------------------- loading
def _expand_features(
    raw: Any, feature_set: str, hierarchy: FeatureHierarchy, ctx: str
) -> frozenset[str]:
    """Expand a ``when`` value (name / list / ``@class``) to a concrete feature set."""
    values = [raw] if isinstance(raw, str) else list(raw)
    out: set[str] = set()
    for value in values:
        if value.startswith("@"):
            node = value[1:]
            if node not in hierarchy.features(feature_set):
                raise SpecError(
                    f"{ctx}: class @{node} is not a feature in feature set {feature_set!r}"
                )
            out |= hierarchy.descendants(feature_set, node, include_root=True)
        else:
            if value not in hierarchy.features(feature_set):
                raise SpecError(f"{ctx}: feature {value!r} is not in feature set {feature_set!r}")
            out.add(value)
    return frozenset(out)


def _check_when_order(keys: Sequence[str], precedence: tuple[str, ...], ctx: str) -> None:
    """Require ``when`` keys to appear in precedence order (readability lint, M2)."""
    positions = [precedence.index(k) for k in keys]  # KeyError-equivalent handled by caller
    if positions != sorted(positions):
        raise SpecError(
            f"{ctx}: 'when' keys {list(keys)} must be written in precedence order {list(precedence)}"
        )


def _parse_emit(raw: Any, precedence: tuple[str, ...], ctx: str) -> Emit:
    if isinstance(raw, str):
        if raw == COMPOSITE:
            return Emit("composite", of="when", sep=_DEFAULT_SEP["when"])
        if raw not in precedence:
            raise SpecError(
                f"{ctx}: emit featureset {raw!r} is not in precedence {list(precedence)}"
            )
        return Emit("featureset", value=raw)
    if "literal" in raw:
        return Emit("literal", value=str(raw["literal"]))
    composite = raw["composite"]
    if isinstance(composite, list):
        featuresets = tuple(composite)
        for fs in featuresets:
            if fs not in precedence:
                raise SpecError(
                    f"{ctx}: composite featureset {fs!r} is not in precedence {list(precedence)}"
                )
        return Emit(
            "composite",
            of="explicit",
            featuresets=featuresets,
            sep=raw.get("sep", _DEFAULT_SEP["explicit"]),
        )
    return Emit("composite", of=composite, sep=raw.get("sep", _DEFAULT_SEP[composite]))


def load_spec(data: Mapping[str, Any], hierarchy: FeatureHierarchy) -> ResolutionSpec:
    """Validate ``data`` (structurally + against ``hierarchy``) and build a spec."""
    try:
        jsonschema.validate(data, SPEC_SCHEMA)
    except jsonschema.ValidationError as exc:
        raise SpecError(f"spec is not well-formed: {exc.message}") from exc

    name = data["name"]
    precedence = tuple(data["precedence"])
    if len(set(precedence)) != len(precedence):
        raise SpecError(f"spec {name!r}: duplicate featureset in precedence {list(precedence)}")
    known = hierarchy.feature_sets()
    for fs in precedence:
        if fs not in known:
            raise SpecError(f"spec {name!r}: unknown feature set {fs!r} in precedence")

    rules: list[Rule] = []
    for i, rule in enumerate(data.get("rules", [])):
        ctx = f"spec {name!r} rule {i}"
        when_raw: Mapping[str, Any] = rule.get("when", {})
        for fs in when_raw:
            if fs not in precedence:
                raise SpecError(f"{ctx}: 'when' feature set {fs!r} is not in precedence")
        _check_when_order(list(when_raw), precedence, ctx)
        when = {fs: _expand_features(val, fs, hierarchy, ctx) for fs, val in when_raw.items()}
        rules.append(Rule(when=when, emit=_parse_emit(rule["emit"], precedence, ctx)))

    return ResolutionSpec(name=name, precedence=precedence, rules=tuple(rules))


def load_spec_file(path: str | Path, hierarchy: FeatureHierarchy) -> ResolutionSpec:
    """Load a YAML resolution spec from ``path``, validated against ``hierarchy``."""
    import yaml

    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, Mapping):
        raise SpecError(f"{path}: top-level YAML must be a mapping")
    return load_spec(data, hierarchy)


#: Built-in presets ship under the ``karyoscope_analysis.presets`` data package.
_PRESET_PACKAGE = "karyoscope_analysis.presets"


def builtin_preset_names() -> list[str]:
    """Names of the built-in resolution presets (the legacy merge_beds priority modes)."""
    return sorted(
        p.name[: -len(".yaml")]
        for p in files(_PRESET_PACKAGE).iterdir()
        if p.name.endswith(".yaml")
    )


def load_builtin_preset(name: str, hierarchy: FeatureHierarchy) -> ResolutionSpec:
    """Load a built-in preset by name, validated against ``hierarchy``."""
    import yaml

    resource = files(_PRESET_PACKAGE) / f"{name}.yaml"
    if not resource.is_file():
        raise SpecError(f"unknown preset {name!r}; available: {builtin_preset_names()}")
    return load_spec(yaml.safe_load(resource.read_text()), hierarchy)
