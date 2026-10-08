# Agent instructions

## Task tracking

- Read [TASK_RULES.md](TASK_RULES.md) when assessing or resuming a task. Automatically create or reuse a canonical record for substantial tasks according to its recovery criteria, and maintain it at meaningful checkpoints.
- Read the task index and matching record when present. Use an explicitly supplied canonical record instead of duplicating it; keep repository-owned tracking suitable for public independent use.
- Clarify consequential scope or visibility uncertainty while continuing independent work. Routine record creation does not require separate permission.

## Specifications

- Before changing behavior, read this repository's `specifications.md` if it exists and the relevant documentation. If it is missing, use existing code and documentation and report consequential uncertainty; do not invent requirements.
- Keep requirements and architectural contracts in the specification, and working instructions in this file.
- Update the specification when an authorized change affects documented behavior, public interfaces, compatibility guarantees, or significant design decisions.
- Keep documentation and development instructions self-contained for an independent checkout. Do not reference private consuming repositories or require access to their specifications or infrastructure.

## Changes and validation

- Preserve unrelated local changes and generated artifacts. Follow existing conventions and keep changes scoped to the task.
- Consider model-format compatibility, public APIs and ABI where applicable, expression evaluation, derivatives, suffixes, and solution output when changing library behavior.
- Use this repository's build and test configuration. Validate affected library variants and relevant numerical, compatibility, and regression behavior; a consuming project's build does not replace component validation.
- Add regression coverage for substantive correctness changes when practical. Use numerical tolerances appropriate to the behavior being checked.
- Report what changed, what was checked, and any validation limitations.
- Keep credentials, confidential models, and private consumer details out of documentation and test fixtures.
