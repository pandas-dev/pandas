///
layout: pandas-dev
title: pandas
author: vgauraha62
///

import pandas as pd
import pytest
import logging

# Starting with logging configuration, to get the insights on the console, important!!
logging.basicConfig(level=logging.INFO)

def test_timestamp_hash_equality_on_dst():
    logging.info("Testing timestamp hash equality on Daylight Saving Time (DST) transition")
    
    # This is a DDaylight Saving Time (DST) transition for the America/Los_Angeles time zone
    dt_str = "2023-11-05 01:00-08:00"
    tz_str = "America/Los_Angeles"

    ts1 = pd.Timestamp(dt_str, tz=tz_str)
    logging.info(f"Created timestamp ts1: {ts1}")
    
    ts2 = ts1 + pd.Timedelta(hours=0)  # This here, should create the same timestamp
    logging.info(f"Created timestamp ts2: {ts2}")

    # Now we verify that the timestamps compare equal
    assert ts1 == ts2, "Timestamps are not considered equal when they should be."
    logging.info("Timestamps are equal")

    # Important to Convert to UTC before comparing hash values to avoid Daylight Saving Time (DST) issues
    ts1_utc = ts1.tz_convert('UTC')
    logging.info(f"Converted ts1 to UTC: {ts1_utc}")
    
    ts2_utc = ts2.tz_convert('UTC')
    logging.info(f"Converted ts2 to UTC: {ts2_utc}")

    # Finally we erify that their hash values are the same, to successfully achieve testing
    assert hash(ts1_utc) == hash(ts2_utc), "Hashes of equal Timestamps do not match after normalization."
    logging.info("Hashes of timestamps are equal")