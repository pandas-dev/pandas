# Pandas to_offset fix
def to_offset(freq):
    """
    Convert frequency string/object to offset.
    Raises ValueError for non-string/non-Offset input.
    """
    if isinstance(freq, str):
        # Process string frequency
        pass
    elif hasattr(freq, 'delta'):
        # Process Offset object
        return freq
    else:
        # Raise ValueError (not TypeError) for invalid input
        raise ValueError(f"Frequency must be string or Offset, got {type(freq).__name__}")
