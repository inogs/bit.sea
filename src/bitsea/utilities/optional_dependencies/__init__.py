class OptionalDependencyMissing(UserWarning):
    pass


# In principle, OptionalDependencyMissing could also have been an
# "ImportWarning", but since they are ignored by default, we use the
# UserWarning
