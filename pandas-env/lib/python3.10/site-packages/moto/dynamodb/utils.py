from typing import List


def extract_duplicates(input_list: List[str]) -> List[str]:
    return [item for item in input_list if input_list.count(item) > 1]
