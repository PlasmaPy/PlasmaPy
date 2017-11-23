"""
plasmapy.errors
===============

Custom Errors and Warning names to improve readability
"""

class PhysicsError(Exception):
    """Error for using a worrisome physics value"""
    pass

class PhysicsWarning(Warning):
    """Warning for using a mildly worrisome physics value"""
    pass

class RelativityWarning(PhysicsWarning):
    """Warning for approaching or exceeding the speed of light"""
    pass