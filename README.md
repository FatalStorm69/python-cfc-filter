# Python CFC Filter

A small Python implementation of an SAE J211-1:1995 Appendix C-style digital Channel Frequency Class (CFC) filter. It calculates the second-order coefficients and uses SciPy's `filtfilt` to apply forward-backward (zero-phase) filtering.

## Supported classes

CFC 60, 180, 600, and 1000. The same Appendix C coefficient equations are used for all four classes for consistency with a MATLAB-style implementation. **SAE J211-1:1995 specifically recommends Appendix C for CFC 60/180; using this method for CFC 600/1000 does not by itself establish conformance.**

## Install

```bash
python -m pip install -r requirements.txt
```

## Usage

```python
from cfc_filter import CFCFilter

filter_ = CFCFilter(cfc_class=60, sampling_frequency=10000)
filtered = filter_.filter_signal(raw_acceleration)  # 1D samples in time order
filter_.verify_frequency_response()                  # prints reference-point checks
```

`raw_acceleration` should be a one-dimensional array-like sequence of numeric samples. Use the **actual acquisition sampling rate**, not an arbitrary rate chosen to obtain a passing response.

Run `python example.py` for a complete synthetic-signal example.

## Frequency-response verification

`verify_frequency_response()` evaluates the **theoretical digital two-pass** attenuation at the reference frequencies using `scipy.signal.freqz`. It prints PASS/FAIL for each evaluable point. Frequencies at or above Nyquist (`sampling_frequency / 2`) are labeled `NOT EVALUABLE`, not PASS or FAIL.

In development tests at 10 kHz, CFC 60, 180 and 600 passed the reference-point checks. CFC 1000 failed the FG lower-bound check at 3496 Hz, and FJ (6442 Hz) could not be evaluated. At 20 kHz, the CFC 1000 reference points passed. These are checks of this particular digital implementation, **not proof of SAE compliance of a complete acquisition or measurement system**.

## Notes

- The coefficients follow the SAE J211-1:1995 Appendix C equations with the `2.0775` frequency factor.
- SciPy's `filtfilt` uses its own endpoint treatment; results near signal boundaries may differ from other implementations.
- Input data should be sampled uniformly and contain enough samples for `filtfilt`.
- This project is not an official SAE implementation, certification, or validation of a measurement system.
- SAE standards are copyrighted; obtain the original standard from an authorized source.
