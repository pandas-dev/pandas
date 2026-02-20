"""Measure document reading durations."""

from __future__ import annotations

import json
import time
from itertools import islice
from operator import itemgetter
from types import NoneType
from typing import TYPE_CHECKING

import sphinx
from sphinx.domains import Domain
from sphinx.locale import __
from sphinx.util import logging

if TYPE_CHECKING:
    from collections.abc import Collection, Set
    from pathlib import Path
    from typing import TypedDict

    from docutils import nodes

    from sphinx.application import Sphinx

    class _DurationDomainData(TypedDict):
        reading_durations: dict[str, float]


logger = logging.getLogger(__name__)


class DurationDomain(Domain):
    """A domain for durations of Sphinx processing."""

    name = 'duration'

    @property
    def reading_durations(self) -> dict[str, float]:
        return self.data.setdefault('reading_durations', {})

    def note_reading_duration(self, duration: float) -> None:
        self.reading_durations[self.env.current_document.docname] = duration

    def warn_reading_duration(self, duration: float, duration_limit: float) -> None:
        logger.warning(
            __('Reading duration %.3fs exceeded the duration limit %.3fs'),
            duration,
            duration_limit,
            type='duration',
            location=self.env.docname,
        )

    def clear(self) -> None:
        self.reading_durations.clear()

    def clear_doc(self, docname: str) -> None:
        self.reading_durations.pop(docname, None)

    def merge_domaindata(  # type: ignore[override]
        self, docnames: Set[str], otherdata: _DurationDomainData
    ) -> None:
        other_reading_durations = otherdata.get('reading_durations', {})
        docnames_set = frozenset(docnames)
        for docname, duration in other_reading_durations.items():
            if docname in docnames_set:
                self.reading_durations[docname] = duration


def on_builder_inited(app: Sphinx) -> None:
    """Initialize DurationDomain on bootstrap.

    This clears the results of the last build.
    """
    domain = app.env.domains['duration']
    domain.clear()


def on_source_read(app: Sphinx, docname: str, content: list[str]) -> None:
    """Start to measure reading duration."""
    app.env.current_document.reading_started_at = time.monotonic()


def on_doctree_read(app: Sphinx, doctree: nodes.document) -> None:
    """Record a reading duration."""
    duration = time.monotonic() - app.env.current_document.reading_started_at
    domain = app.env.domains['duration']
    domain.note_reading_duration(duration)

    duration_limit: float | None = app.config.duration_limit
    if duration_limit is not None and duration > duration_limit:
        domain.warn_reading_duration(duration, duration_limit)


def on_build_finished(app: Sphinx, error: Exception) -> None:
    """Display duration ranking on the current build."""
    domain = app.env.domains['duration']
    if not domain.reading_durations:
        return

    # Get default options and update with user-specified values
    if app.config.duration_print_total:
        _print_total_duration(domain.reading_durations.values())

    if app.config.duration_print_slowest:
        _print_slowest_durations(
            domain.reading_durations, app.config.duration_n_slowest
        )

    if write_json := app.config.duration_write_json:
        _write_json_durations(domain.reading_durations, app.outdir / write_json)


def _print_total_duration(durations: Collection[float]) -> None:
    logger.info('')
    logger.info(
        __('====================== total reading duration ==========================')
    )

    n_files = len(durations)
    s = 's' if n_files != 1 else ''
    minutes, seconds = divmod(sum(durations), 60)
    logger.info(
        __('Total time reading %d file%s: %dm %.3fs'), n_files, s, minutes, seconds
    )


def _print_slowest_durations(durations: dict[str, float], n_slowest: int) -> None:
    sorted_durations = sorted(durations.items(), key=itemgetter(1), reverse=True)
    n_slowest = n_slowest or len(sorted_durations)
    n_slowest = min(n_slowest, len(sorted_durations))

    logger.info('')
    logger.info('')
    logger.info(
        __('====================== slowest reading durations =======================')
    )
    for docname, duration in islice(sorted_durations, n_slowest):
        logger.info(__('%.3fs %s'), duration, docname)

    logger.info('')


def _write_json_durations(durations: dict[str, float], out_file: Path) -> None:
    durations = {k: round(v, 3) for k, v in durations.items()}
    out_file.parent.mkdir(parents=True, exist_ok=True)
    durations_json = json.dumps(durations, ensure_ascii=False, indent=4, sort_keys=True)
    out_file.write_text(durations_json, encoding='utf-8')


def setup(app: Sphinx) -> dict[str, bool | str]:
    app.add_domain(DurationDomain)
    app.connect('builder-inited', on_builder_inited)
    app.connect('source-read', on_source_read)
    app.connect('doctree-read', on_doctree_read)
    app.connect('build-finished', on_build_finished)

    app.add_config_value('duration_print_total', True, '', types=frozenset({bool}))
    app.add_config_value('duration_print_slowest', True, '', types=frozenset({bool}))
    app.add_config_value('duration_n_slowest', 5, '', types=frozenset({int}))
    app.add_config_value(
        'duration_write_json',
        'sphinx-reading-durations.json',
        '',
        types=frozenset({str, NoneType}),
    )
    app.add_config_value(
        'duration_limit', None, '', types=frozenset({float, int, NoneType})
    )

    return {
        'version': sphinx.__display_version__,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
