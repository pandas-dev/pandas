from typing import Any, Dict, Optional


class BadSegmentException(Exception):
    def __init__(
        self,
        seg_id: Optional[str] = None,
        code: Optional[str] = None,
        message: Optional[str] = None,
    ):
        self.id = seg_id
        self.code = code
        self.message = message

    def __repr__(self) -> str:
        return f"<BadSegment {self.id}-{self.code}-{self.message}>"

    def to_dict(self) -> Dict[str, Any]:
        result: Dict[str, Any] = {}
        if self.id is not None:
            result["Id"] = self.id
        if self.code is not None:
            result["ErrorCode"] = self.code
        if self.message is not None:
            result["Message"] = self.message

        return result
