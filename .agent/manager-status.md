# Agent Manager Status

## Current Focus
Foreground-marginalized CMB bandpowers integration for `lens_only` mode.

## Session History

### 2026-01-06: Initial setup
- Initialized agent management files
- Reviewed codebase structure and existing `lens_only` implementation
- Identified key code paths and existing test coverage

### 2026-01-06: Foreground marginalization focus
- User clarified: focus on foreground marginalized step
- Source data: `/home/jiaqu/DR6-ACT-lite/` contains DR6 foreground-marginalized CMB bandpowers
- Task: Read README, extract bandpowers and fiducial data as needed

## Outstanding Tasks
1. **[IN PROGRESS]** Explore `/home/jiaqu/DR6-ACT-lite/` — read README, identify bandpower files and fiducial data
2. Extract/integrate foreground-marginalized CMB bandpowers into dr6plus_lenslike
3. Add/improve tests for `lens_only=True` pathway with foreground marginalization
4. Ensure backward compatibility

## Notes
- Existing tests in `tests/test_dr6plus_lenslike.py` use mocked data
- `lens_only=True` currently loads CMB-marginalized covariance and skips likelihood corrections
- Previous session died when accessing `/home/jiaqu/DR6-ACT-lite/` — proceed carefully
- Foreground marginalization is the key focus for the CMB 2pt marginalization approach
