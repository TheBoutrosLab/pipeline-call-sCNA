# Description
<!-- Briefly describe the changes included in this pull request and provide paths
to the test cases below. Start with "Closes #..." when appropriate. -->

### Closes #...

## Testing Results
<!-- Delete unused cases. Include NFTest results when tests are available. -->

- CNV_FACETS
    - sample:     <!-- e.g. paired A-mini sample -->
    - input YAML: <!-- path/to/input-call-sCNA.yaml -->
    - config:     <!-- path/to/sample.config -->
    - output:     <!-- path/to/CNV_FACETS output -->
- Battenberg
    - sample:     <!-- e.g. paired A-mini WGS sample -->
    - input YAML: <!-- path/to/input-call-sCNA.yaml -->
    - config:     <!-- path/to/sample.config -->
    - output:     <!-- path/to/Battenberg output -->
- NFTest
    - output: <!-- path/to/output -->
    - log:    <!-- path/to/log -->
    - cases:  <!-- list of test cases run -->

# Checklist
<!-- Confirm each applicable item by replacing [ ] with [x]. -->

- [ ] I have read the [code review guidelines](https://solid-adventure-l491og6.pages.github.io/latest/code-review-guidelines) and the [code review best-practices checklist](https://solid-adventure-l491og6.pages.github.io/latest/code-review-best-practices).

- [ ] I have reviewed the [Nextflow pipeline standards](https://solid-adventure-l1qkwg6.pages.github.io/latest/nextflow-pipeline-standardization).

- [ ] The branch name follows the repository standards: `[AD username]-[brief description]`.

- [ ] I have set up or verified the branch protection rules before opening this pull request.

- [ ] I have added my name to the `manifest` contributors in `nextflow.config`, am already listed, or do not wish to be listed. This acknowledgement is optional.

- [ ] I have documented the changes in `CHANGELOG.md` under the next release or unreleased section and updated the date.

- [ ] I have updated the `manifest` version in `nextflow.config` following [semantic versioning](https://semver.org/), or the version has already been updated. Leave this unchecked and discuss it in the pull request if unsure.

- [ ] I have tested each affected caller on an appropriate paired sample and recorded the test inputs and outputs above.
