"""SAE J211-style CFC digital filtering and frequency-response checks.

Uses Appendix C coefficient equations and SciPy forward-backward filtering.
This module does not certify an acquisition system's SAE compliance.
"""

from typing import ClassVar
import math

import numpy as np
from scipy.signal import filtfilt, freqz


class CFCFilter:
    """Two-pass, zero-phase digital filter for CFC 60, 180, 600 or 1000."""

    # Reference frequencies (Hz) for response verification.
    CFC_REQUIREMENTS: ClassVar[dict] = {
        60: {"FL": 0.1, "FH": 60, "FC": 75, "FN": 100, "FD": 130, "FE": 160, "FG": 452},
        180: {"FL": 0.1, "FH": 180, "FC": 225, "FN": 300, "FD": 390, "FE": 480, "FG": 1310},
        600: {"FL": 0.1, "FH": 600, "FN": 1000, "2*FH": 1200, "FG": 2119, "FJ": 3865},
        1000: {"FL": 0.1, "FH": 1000, "FN": 1600, "2*FH": 2000, "FG": 3496, "FJ": 6442},
    }

    # Lower and upper attenuation limits (dB); None means no bound.
    CFC_LIMITS: ClassVar[dict] = {
        60: {
            "FL": {"upper": 0.5, "lower": -0.5},
            "FH": {"upper": 0.5, "lower": -1.0},
            "FC": {"upper": -0.3, "lower": -1.8},
            "FN": {"upper": -1.8, "lower": -3.8},
            "FD": {"upper": -5.2, "lower": -8.2},
            "FE": {"upper": -9.2, "lower": -13.2},
            "FG": {"upper": -40.0, "lower": -48.3},
        },
        180: {
            "FL": {"upper": 0.5, "lower": -0.5},
            "FH": {"upper": 0.5, "lower": -1.0},
            "FC": {"upper": -0.3, "lower": -1.8},
            "FN": {"upper": -1.8, "lower": -3.8},
            "FD": {"upper": -5.2, "lower": -8.2},
            "FE": {"upper": -9.2, "lower": -13.2},
            "FG": {"upper": -40.0, "lower": -48.3},
        },
        600: {
            "FL": {"upper": 0.5, "lower": -0.5},
            "FH": {"upper": 0.5, "lower": -1.0},
            "FN": {"upper": 0.5, "lower": -4.0},
            "2*FH": {"upper": 0.5, "lower": None},
            "FG": {"upper": None, "lower": -30.0},
            "FJ": {"upper": -40.0, "lower": float("-inf")},
        },
        1000: {
            "FL": {"upper": 0.5, "lower": -0.5},
            "FH": {"upper": 0.5, "lower": -1.0},
            "FN": {"upper": 0.5, "lower": -4.0},
            "2*FH": {"upper": 0.5, "lower": None},
            "FG": {"upper": None, "lower": -30.0},
            "FJ": {"upper": -40.0, "lower": float("-inf")},
        },
    }

    def __init__(self, cfc_class: int, sampling_frequency: float):
        if cfc_class not in self.CFC_REQUIREMENTS:
            raise ValueError("Unsupported CFC class")
        if not math.isfinite(sampling_frequency) or sampling_frequency <= 0:
            raise ValueError("Sampling frequency must be positive and finite")
        self.cfc_class = cfc_class
        self.sampling_frequency = sampling_frequency

    def get_frequencies(self):
        """Return reference frequencies for the selected CFC class."""
        return self.CFC_REQUIREMENTS[self.cfc_class]

    def get_limits(self):
        """Return reference attenuation limits for the selected class."""
        return self.CFC_LIMITS[self.cfc_class]

    def calculate_coefficients(self):
        """Return numerator B and denominator A for a second-order filter."""
        dt = 1 / self.sampling_frequency
        wd = 2 * math.pi * self.cfc_class * 2.0775
        wa = math.sin(wd * dt / 2) / math.cos(wd * dt / 2)

        # Appendix C coefficient equations (single-pass section).
        denominator = 1 + math.sqrt(2) * wa + wa**2
        a0 = wa**2 / denominator
        a1 = 2 * a0
        a2 = a0
        b1 = -2 * (wa**2 - 1) / denominator
        b2 = (-1 + math.sqrt(2) * wa - wa**2) / denominator

        # scipy.signal uses negative feedback coefficients in A.
        B = [a0, a1, a2]
        A = [1, -b1, -b2]
        return B, A

    def filter_signal(self, signal):
        """Filter a 1D signal forward and backward for zero phase."""
        B, A = self.calculate_coefficients()
        return filtfilt(B, A, signal)

    def verify_frequency_response(self):
        """Print digital two-pass attenuation checks at reference frequencies.

        Points at or above Nyquist are reported as NOT EVALUABLE.
        This checks the digital filter only, not the measurement chain.
        """
        B, A = self.calculate_coefficients()
        frequencies = self.get_frequencies()
        limits = self.get_limits()
        nyquist = self.sampling_frequency / 2

        for label, frequency in frequencies.items():
            if frequency >= nyquist:
                print(f"{label:5} {frequency:8.1f} Hz  NOT EVALUABLE (Nyquist: {nyquist:g} Hz)")
                continue

            _, h = freqz(B, A, worN=np.array([frequency]), fs=self.sampling_frequency)
            # Two passes square the magnitude, doubling attenuation in dB.
            attenuation = float(40 * np.log10(np.abs(h[0])))
            lower = limits[label]["lower"]
            upper = limits[label]["upper"]
            lower_ok = lower is None or attenuation >= lower
            upper_ok = upper is None or attenuation <= upper
            result = "PASS" if lower_ok and upper_ok else "FAIL"
            print(f"{label:5} {frequency:8.1f} Hz  {attenuation:9.3f} dB  {result}")
