
"""
PAUdm constants definition

@author: anne bauer, santi serrano

Defines constant such as flag bit meanings.

"""

# Mask flags
flag_crowded = 1
flag_merged = 2
flag_halo = 4
flag_truncated = 8
flag_deblended = 16
flag_crosstalk = 32
flag_scatterlight = 64
flag_extinction = 128
flag_photo_zp = 256
flag_cosmetics = 512  # (image mask)
flag_saturated = 1024  # (image mask)
flag_cosmics = 2048  # (image mask)
flag_vignetted = 4096  # (image mask)
flag_discordant = 8192
flag_edge = 16384
flag_distortion = 32768
flag_noise = 65536
flag_astrometry = 131072


# Operation - Job states
OPEN_LAUNCHED_STATES = ('LAUNCHED', 'SUBMITTED', 'WAITING', 'READY', 'SCHEDULED', 'RUNNING')
OPEN_FAILED_STATES = ('ABORTED', 'CANCELLED', 'LAUNCH_FAIL', 'CODE_FAIL', 'CLEAR_FAIL')
OPEN_STATES =  ('NEW',) + OPEN_LAUNCHED_STATES + OPEN_FAILED_STATES
CLOSED_FAILED_STATES = ('ABORTED_CLOSED', 'CANCELLED_CLOSED', 'LAUNCH_FAIL_CLOSED', 'CODE_FAIL_CLOSED', 'CLEAR_FAIL_CLOSED', 'PARENT_FAIL_CLOSED')
CLOSED_STATES = ('CLEARED_CLOSED',) + CLOSED_FAILED_STATES


# Photometry constants
mag_zp = 26.0

# Match radius for global_objects
match_radius = 1.0/3600.0  # degrees

