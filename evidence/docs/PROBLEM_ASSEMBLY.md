# Problem Assembly

`assemble_problem_spec(...)` does one of two things:

1. returns an assembled `nomodescent.ProblemSpec`
2. returns an `underdetermined_evidence` result with the unresolved observer
   groups and required human decisions

It never silently fills in missing structure.

Assembly currently requires:

- enough observer hypotheses to resolve one observer per protocol
- at least one surviving visible matrix observation
- no conflicting authoritative visible observations for the same observer and
  matrix kind

Encoded observer hypotheses can be promoted into a `ProblemSpec`, but only after
they are either:

- uniquely determined, or
- explicitly selected by the caller
