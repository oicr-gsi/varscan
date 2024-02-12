## 2.4.1 - 2024-02-12
- Added minMemory parameter to tasks using RAM scaling (to handle small chromosomes better)
## 2.4.0 - 2023-12-16
- Added scaling RAM assignment by chromosome for scattered tasks
## 2.3.0 - 2023-06-21
- Assembly-specific settings go into wdl
## 2.2.5 - 2023-06-08
- Increment version to install varscan_matched in Vidarr
## Unreleased
- updated README, explained passing interval file better
## 2.2.4 - 2022-05-25
- Change most output types to non-optional. Classified as a bug fix. Filtered outputs are still optional!
## 2.2.3 - 2022-02-02
- Add timeout and modules param to interval-processing task so that we can use modularized data
## 2.2.2 - 2021-06-02
- This version increment was to resolve installation issues in vidarr
## Unreleased - 2021-11-10
[GP-2891](https://jira.oicr.on.ca/browse/GP-2891) Making RT tests more robust
## 2.2.1 - 2021-02-01
- Increment version to avoid overlap with a compromized installation
## 2.2   - 2021-01-15
- Re-implemented to cut on merging pileups and scatter variant-calling tasks
## 2.1.1 - 2020-04-06
- Added timeout setting to mpileup concatenation step. 5 hr is not enough for some data
## 2.1 - 2020-04-02
- Introduced compression of mpileup data to conserve disk space
## 2.0 - 2020-01-31
- Converting to wdl workflow, adding more output types to mesh well with sequenza
## 1.0 - 2017-08-16
- Initial implementation as a stand-alone workflow
