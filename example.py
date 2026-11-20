"""Minimal example using a synthetic signal (no external test data)."""

import numpy as np
from cfc_filter import CFCFilter

sampling_frequency = 10000
seconds = np.arange(0, 0.2, 1 / sampling_frequency)
raw_acceleration = np.sin(2 * np.pi * 50 * seconds) + 0.2 * np.sin(2 * np.pi * 1500 * seconds)

cfc_filter = CFCFilter(cfc_class=60, sampling_frequency=sampling_frequency)
filtered_acceleration = cfc_filter.filter_signal(raw_acceleration)
print(f"Raw peak: {raw_acceleration.max():.4f}")
print(f"Filtered peak: {filtered_acceleration.max():.4f}")
cfc_filter.verify_frequency_response()
