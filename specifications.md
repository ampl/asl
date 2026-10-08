# ASL specifications

Documentation baseline: 2026-10-08. Initial specification of the checked-out public AMPL Solver Library. **Current** statements describe inspected implementation/documentation; **proposed** statements describe acceptance or maintenance policy for review. Open decisions are in section 8. Build configurations and historical documentation are not blanket support or ABI guarantees.

## 1. Project overview

ASL is a public collection of libraries and supporting tools that help solver interfaces read AMPL model files, evaluate expressions and derivatives, exchange solver options/suffixes, and write solution results.

Consumers are solver-interface developers, applications using ASL's model/evaluation APIs, and users of associated examples and solver interfaces. Expected value is reusable model access and evaluation infrastructure that can be built and used independently.

Scope includes the `asl`/`asl2` source families, build variants, optional C++ wrappers, function-linking and Fortran interoperability support, examples, and distributed interface source. Optimization algorithms and consumer-specific packaging remain outside the core library contract.

## 2. Requirements

### 2.1 Functional requirements

**Current contracts and capabilities:**

- Read supported AMPL NL model data and expose problem structure, variables, bounds, objectives, constraints, and applicable suffix information through the chosen library/API family.
- Evaluate model functions and supported derivatives, including objective gradients, constraint Jacobians, and Hessian-related operations as exposed by the selected API and read mode.
- Provide solver-facing option, function-linking, suffix, and SOL output facilities. Associated examples demonstrate integration and evaluation.
- Build `asl` and `asl2`, with optional shared, runtime-linkage, multithread, C++ wrapper, and supporting/example variants selected by CMake.
- Provide optional f2c support and documented calling conventions for relevant Fortran interfaces.

**Proposed acceptance contracts:**

- Supported input is interpreted consistently with its documented format/API conventions. Invalid input and evaluation-domain failures must not be mistaken for successful evaluation.
- Values and derivatives agree with analytic or independently verified references within explicitly chosen tolerances, including sparse structure and indexing where relevant.
- Solution writing preserves the intended association of values/statuses/suffixes with model entities. Output availability depends on the calling solver and API; ASL does not itself guarantee optimizer convergence.

The human manual is section 2.8, with links to API/example documentation.

### 2.2 Data requirements

Primary input is NL model data, together with applicable auxiliary files, user-defined functions, solver options, and suffixes. Output includes in-memory model/evaluation data, SOL results, and example-specific files such as `gjh` derivative output.

Generated platform headers describe arithmetic and system characteristics. Their configuration must match the target architecture; they are build products, not portable model data. Cross-compilation requires attention to the target's arithmetic assumptions.

Validation fixtures should cover supported format/read modes, variable/constraint ordering, sparsity, derivative values, and failures. Use public synthetic or publishable models. No freshness requirement or persistent user-data store applies to the core library.

### 2.3 Non-functional requirements and compatibility

**Current:** CMake defines static/shared outputs, optional Windows dynamic-runtime libraries, multithread variants, and optional C++ wrappers. The checked-in CI declares Linux, MSVC, MinGW, and macOS builds with selected variants; historical runner labels do not establish current availability or complete support.

**Proposed:** maintain documented format/API behavior and validate numerical changes across relevant families/read modes/platforms. Thread-safety is specific to APIs, contexts, variants, and caller behavior; an OpenMP-enabled build alone does not promise unrestricted sharing of model state. No blanket thread-safety, API/ABI stability, latency, or supported platform matrix is established by this initial specification.

### 2.4 Operational and quality requirements

| Validation layer | Source | Scope and limitation |
| --- | --- | --- |
| Standalone build/variant checks | [CMakeLists.txt](CMakeLists.txt), [.github/workflows/cmake.yml](.github/workflows/cmake.yml) | CI builds selected platforms/variants; inspected workflow contains no numerical test execution. |
| Evaluation and derivative exercises | [src/examples](src/examples), including `evalchk` and `gjh` | Useful executable entry points; building examples is not proof that numerical cases ran. |
| Caller/integration regressions | Public applications and solver-interface fixtures | Validate real model reads/evaluation/result output with explicit consumer settings; do not replace independent library checks. |
| Proposed numerical regression suite | To be selected/organized by maintainers | Analytic/reference values, gradients, Jacobians, Hessian-related operations, sparsity/indexing, invalid inputs, and relevant variants. |

No top-level CTest suite is registered in the inspected ASL CMake configuration. **Proposed:** record numerical tests separately from build checks and add reproducible cases for correctness fixes. Finite-difference checks can supplement analytic references on well-conditioned, smooth cases; discontinuities and domain boundaries need deliberate expectations and tolerances.

Logs and example output should diagnose read/evaluation failures without unnecessarily exposing confidential model contents. The core library is not an always-on service; availability/monitoring/backup obligations belong to applications that embed it.

### 2.5 External integrations

| Integration | Purpose and constraints |
| --- | --- |
| C/C++ toolchains and platform libraries | Compile/link selected libraries and wrappers; arithmetic and runtime choices must match the target. |
| OpenMP/toolchain support | Used by multithread variants; availability and runtime distribution depend on platform/toolchain. |
| f2c and related support | Optional Fortran calling/support path and examples; callers must use compatible conventions. |
| User-defined function libraries | Optional externally supplied evaluation functions; API and lifetime/error behavior must match the calling ASL context. |
| AMPL or another compatible NL producer | Provides model input and consumes results; not required to compile the core libraries. |
| Vendor solver libraries | Needed only for particular solver interface source, not the default core-library build. |

Public consumers may embed ASL or use installed libraries; independent setup must not depend on private consuming repositories.

### 2.6 Technology constraints

[CMakeLists.txt](CMakeLists.txt) and [support/cmake](support/cmake) define architecture checks, arithmetic-header generation, optimization, variants, and install targets. Historical make/configure workflows are described in the source READMEs. Use current build declarations for exact option spelling and behavior.

`GENERATE_ARITH` invokes host arithmetic detection when enabled; the CMake comments explicitly warn that this path is not for cross-compiling. Static/shared linkage, Windows runtime selection, and installed headers must be compatible with the consumer's build. Dependency versions and generator requirements belong in maintained build configuration.

### 2.7 Explicit non-goals

- Implement a general optimization engine or guarantee a solver's convergence/results.
- Promise that `asl` and `asl2` are interchangeable for every client or read/evaluation mode.
- Guarantee thread safety merely because a multithread variant was selected.
- Define downstream application authentication, deployment, licensing, or user-data persistence.

### 2.8 Human manual: standalone build and usage

These source-derived commands have not been executed as part of preparing this specification. Run from an independent ASL checkout with CMake and an appropriate compiler, using a new build tree:

```console
cmake -S . -B build-spec -DCMAKE_BUILD_TYPE=Release -DBUILD_ASL_EXAMPLES=ON
cmake --build build-spec --config Release
```

Current CMake uses `BUILD_ASL_EXAMPLES`; the README and some CI commands use the historical name `BUILD_EXAMPLES`. Use the current option to request examples. Examples enable f2c support and additional targets. To build only the default libraries, omit the examples option.

Optional configurations include `BUILD_CPP`, `BUILD_SHARED_LIBS`, `BUILD_MT_LIBS`, and, on Windows, `BUILD_DYNRT_LIBS`. Choose only variants compatible with the consumer's API/runtime. A selected configuration may require additional compiler/runtime support. For installation into a local staging directory:

```console
cmake --install build-spec --config Release --prefix build-spec/install
```

The install rules export libraries, generated/public headers for the chosen source families, and CMake target metadata. Outputs normally reside under `lib/` and `bin/`; multi-configuration generators may add configuration subdirectories. See [README.md](README.md) for platform examples, reconciling historical switches with [CMakeLists.txt](CMakeLists.txt).

To use ASL, include the headers for the chosen API family, link its compatible target/library, allocate/read the appropriate model context, perform supported evaluation or solver integration, write applicable results, and release resources according to the API. Detailed calling conventions and examples are in [src/solvers/README](src/solvers/README), [src/solvers2/README](src/solvers2/README), their headers, and [src/examples/README](src/examples/README).

For a concrete derivative-output workflow, build examples, place `gjh` on the executable search path, and follow [README.gjh](src/examples/README.gjh) with a valid AMPL model. It computes gradient/Jacobian/Hessian information and emits a file for inspection; this is evaluation, not optimization. Read the example's output to locate its generated file.

The inspected standalone CMake project does not register CTest numerical tests. Validate changed evaluations using appropriate reference fixtures and example/client programs, and record commands, inputs, expected values, tolerances, and results. Do not report a successful build as a numerical test pass.

## 3. Design choices

### 3.1 Project structure and ownership

| Location | Responsibility |
| --- | --- |
| `src/solvers/` | ASL source family, public headers, NL/SOL handling, and evaluation infrastructure |
| `src/solvers2/` | ASL2 source family and corresponding headers/infrastructure |
| `src/cpp/` | Optional C++ wrappers |
| `src/f2c/`, `src/funclink/` | Optional Fortran support and function-linking facilities |
| `src/examples/` | Evaluation and solver-interface examples and related documentation |
| `solvers/` | Distributed solver-interface source; engine dependencies vary |
| `support/cmake/`, root CMake | Build/architecture configuration and library installation |
| `.github/workflows/` | Standalone CI build definitions |
| Generated `include/`, `bin/`, `lib/` in build trees | Platform headers, executables, and libraries; generated artifacts |

### 3.2 Interfaces and dependency direction

Clients select an ASL API/source family and compatible build variant. They own the solver algorithm and application behavior; ASL owns supported model access/evaluation and related library interfaces. Optional wrappers and examples depend on core libraries. Optional user-function and solver integrations introduce explicit external dependencies rather than making every ASL consumer require every engine.

Documentation and development instructions are self-contained for public independent use. Consuming projects reference ASL contracts and maintain their own integration requirements separately.

### 3.3 Important decisions

- **Separate ASL/ASL2 families (current):** maintained source families expose distinct implementation/API paths. Consequence: validate the family/read modes actually used rather than assuming one family's result proves the other.
- **Target-specific arithmetic generation (current):** generated headers encode platform arithmetic characteristics. Consequence: cross-builds must not substitute host characteristics without validation.
- **Selectable library variants (current):** runtime linkage, shared/static, multithread, and wrapper choices support different clients. Consequence: matching headers, runtime, and concurrency assumptions are part of integration correctness.
- **Examples alongside libraries (current):** example programs make API workflows concrete. Consequence: use them as reproducible validation entry points, while recording whether numerical checks actually ran.

Historical rationale and planned changes beyond the documented build comments need maintainer confirmation.

### 3.4 Sources of truth and related documents

| Information | Source |
| --- | --- |
| API/read/evaluation contracts | Headers and source in [src/solvers](src/solvers) and [src/solvers2](src/solvers2), their READMEs |
| Options, variants, generated headers, install targets | [CMakeLists.txt](CMakeLists.txt), [support/cmake](support/cmake) |
| Platform setup | [README.md](README.md), reconciled with current build definitions |
| Executable/API examples | [src/examples](src/examples), [README.gjh](src/examples/README.gjh), [README.f77](src/solvers/README.f77), [README.suf](src/solvers/README.suf) |
| CI build selections | [.github/workflows/cmake.yml](.github/workflows/cmake.yml) |
| License notices | [LICENSE](LICENSE), [LICENSE.2](LICENSE.2), applicable source notices |
| Agent working instructions | [AGENTS.md](AGENTS.md) |

### 3.5 Change coordination

Update this specification and API/example documentation when supported behavior changes. Preserve calling, indexing, arithmetic, and ownership conventions or document intentional changes. Validate affected source families and variants, and communicate migration needs to public consumers without depending on their private infrastructure. Numerical fixes should include a publishable reproducer and expected results.

## 4. Security and privacy

ASL is public source; its license notices are in the repository. Models and derivative/solution exports can contain confidential data. Native parsing/evaluation and externally loaded user functions are important input boundaries.

**Proposed:** use controlled fixtures, review malformed-input and allocation/error paths when changing them, and exercise memory-safety tooling where practical. User-defined function code must be trusted according to the embedding application's execution policy; compiling the library does not provide a sandbox for that code. Exclude secrets, credentials, and private user models from public source/tests. ASL provides no application identity service or confidential-data store. This specification is not a completed security audit.

## 5. Analytics and success criteria

**Proposed acceptance:** intended standalone targets/variants compile and install; public client examples link; relevant model reads and result writes behave as documented; and reference value/derivative cases satisfy declared tolerances and structure/indexing checks. Build passes and numerical validation passes must be recorded separately.

Track tested families, read modes, platforms, variants, passed/failed cases, numerical discrepancies, and performance changes with reproducible metadata. Establish evaluation-time/memory baselines for performance work before choosing thresholds. Global coverage, latency targets, API/ABI guarantees, and numerical tolerance policies remain open decisions.

## 6. Milestones

Proposed ongoing stages, without calendar commitments or completion claims.

| Stage | Deliverable | Acceptance/dependency |
| --- | --- | --- |
| Review baseline | Confirm API/support boundaries and documentation gaps | Maintainer review |
| Establish numerical validation baseline | Publishable reference fixtures and executable checks | Chosen families/read modes, expected values, and tolerances |
| Implement a library fix or extension | Source, documentation, and reproducer | Affected build variants and reference checks pass |
| Validate and communicate a release | Standalone build/install evidence and compatibility notes | Supported matrix and migration policy agreed |

## 7. Use cases

- **Integrate an optimizer:** select an API family, read NL input, evaluate required functions/derivatives, invoke the engine, and write applicable solution information.
- **Inspect derivatives:** use an example or client to compare function/derivative data with an analytic reference, including sparsity and indexing.
- **Build a platform variant:** choose runtime/concurrency/linkage settings compatible with the client and verify installation and linking.
- **Fix a read/evaluation regression:** reduce to a public fixture, record the family/read mode and expected behavior, and validate the affected paths.

## 8. Open decisions and known gaps

- Define supported platforms/toolchains, family/read-mode coverage, API/ABI stability scope, and deprecation policy.
- Confirm thread-safety contracts for selected variants/APIs and allowable shared-context usage.
- Reconcile README/CI `BUILD_EXAMPLES` references with current `BUILD_ASL_EXAMPLES` behavior.
- Establish a maintained numerical regression entry point; inspected top-level CI builds but does not execute such a suite.
- Choose reference models, derivative tolerances, difficult-domain expectations, and performance baselines.
- Confirm ownership and evidence required for standalone release acceptance; historical workflow matrix entries are not proof of current support.
