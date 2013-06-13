
"""
PAUdm constants definition

@author: anne bauer

Defines constant such as flag bit meanings.

"""

# Mask flags
flag_sex_crowded = 1
flag_sex_merged = 2
flag_sex_saturated = 4
flag_sex_truncated = 8
flag_sex_aperture_bad = 16
flag_sex_isophote_bad = 32
flag_sex_deblend_overflow = 64
flag_sex_extract_overflow = 128
flag_cosmetics_noninterp = 256 # Cosmetic non-interpolated
flag_cosmetics_interp = 512 # Cosmetic interpolated
flag_cosmetics_saturated = 1024 # Saturated pixel
flag_cosmics = 2048 # Cosmic Ray
flag_vignetted = 4096
flag_max = 8191 # Not an actual flag, but the maximum possible flag integer given the above flags.



# Operation - Job states
OPEN_LAUNCHED_STATES = ('LAUNCHED', 'SUBMITTED', 'WAITING', 'READY', 'SCHEDULED', 'RUNNING')
OPEN_FAILED_STATES = ('ABORTED', 'CANCELLED', 'LAUNCH_FAIL', 'CODE_FAIL', 'CLEAR_FAIL')
OPEN_STATES =  ('NEW',) + OPEN_LAUNCHED_STATES + OPEN_FAILED_STATES
CLOSED_FAILED_STATES = ('ABORTED_CLOSED', 'CANCELLED_CLOSED', 'LAUNCH_FAIL_CLOSED', 'CODE_FAIL_CLOSED', 'CLEAR_FAIL_CLOSED', 'PARENT_FAIL_CLOSED')
CLOSED_STATES = ('CLEARED_CLOSED',) + CLOSED_FAILED_STATES


# Photometry constants
mag_zp = 26.0

# Match radius for global_objects
match_radius = 1.0/3600.0 # degrees
