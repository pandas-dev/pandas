###############################################################################
#
# ChartTitle - A class for representing Excel chart titles.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

from typing import Any, Dict, Optional


class ChartTitle:
    """
    A class to represent an Excel chart title.

    This class encapsulates all title related properties and methods for the
    chart title and axis titles.
    """

    def __init__(self) -> None:
        """
        Initialize a ChartTitle instance.
        """
        self.font: Optional[Dict[str, Any]] = None
        self.name: Optional[str] = None
        self.formula: Optional[str] = None
        self.data_id: Optional[int] = None
        self.layout: Optional[Dict[str, Any]] = None
        self.overlay: Optional[bool] = None
        self.hidden: bool = False
        self.line: Optional[Dict[str, Any]] = None
        self.fill: Optional[Dict[str, Any]] = None
        self.pattern: Optional[Dict[str, Any]] = None
        self.gradient: Optional[Dict[str, Any]] = None

    def has_name(self) -> bool:
        """
        Check if the title has a text name set.

        Returns:
            True if name has been set.
        """
        return self.name is not None and self.name != ""

    def has_formula(self) -> bool:
        """
        Check if the title has a formula set.

        Returns:
            True if formula has been set.
        """
        return self.formula is not None

    def has_formatting(self) -> bool:
        """
        Check if the title has any formatting properties set.

        Returns:
            True if the title has line, fill, pattern, or gradient formatting.
        """
        has_line = self.line is not None and self.line.get("defined", False)
        has_fill = self.fill is not None and self.fill.get("defined", False)
        has_pattern = self.pattern
        has_gradient = self.gradient

        return has_line or has_fill or has_pattern or has_gradient

    def get_formatting(self) -> Dict[str, Any]:
        """
        Get a dictionary containing the formatting properties.

        Returns:
            A dictionary with line, fill, pattern, and gradient properties.
        """
        return {
            "line": self.line,
            "fill": self.fill,
            "pattern": self.pattern,
            "gradient": self.gradient,
        }

    def is_hidden(self) -> bool:
        """
        Check if the title is explicitly hidden.

        Returns:
            True if title is hidden.
        """
        return self.hidden

    def __repr__(self) -> str:
        """
        Return a string representation of the ChartTitle.
        """
        return (
            f"ChartTitle(\n"
            f"    name = {self.name!r},\n"
            f"    formula = {self.formula!r},\n"
            f"    hidden = {self.hidden!r},\n"
            f"    font = {self.font!r},\n"
            f"    line = {self.line!r},\n"
            f"    fill = {self.fill!r},\n"
            f"    pattern = {self.pattern!r},\n"
            f"    gradient = {self.gradient!r},\n"
            f"    layout = {self.layout!r},\n"
            f"    overlay = {self.overlay!r},\n"
            f"    has_formatting = {self.has_formatting()!r},\n"
            f")\n"
        )
