#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Unit tests for diphtoscan.logging_config module."""

import logging
import pytest
from diphtoscan.logging_config import (
    setup_logging,
    get_logger,
    log_debug,
    log_info,
    log_warning,
    log_error
)


class TestSetupLogging:
    """Tests for setup_logging function."""

    def test_default_logging_level_is_info(self):
        """Default logging should be INFO level."""
        logger = setup_logging()
        assert logger.level == logging.INFO

    def test_verbose_mode_sets_debug_level(self):
        """Verbose mode should set DEBUG level."""
        logger = setup_logging(verbose=True)
        assert logger.level == logging.DEBUG

    def test_quiet_mode_sets_error_level(self):
        """Quiet mode should set ERROR level."""
        logger = setup_logging(quiet=True)
        assert logger.level == logging.ERROR

    def test_logger_has_handler(self):
        """Logger should have exactly one handler after setup."""
        logger = setup_logging()
        assert len(logger.handlers) == 1

    def test_logger_name_is_diphtoscan(self):
        """Logger name should be 'diphtoscan'."""
        logger = setup_logging()
        assert logger.name == 'diphtoscan'

    def test_multiple_setup_calls_clear_handlers(self):
        """Multiple setup calls should not accumulate handlers."""
        setup_logging()
        setup_logging()
        setup_logging()
        logger = get_logger()
        assert len(logger.handlers) == 1

    def test_quiet_overrides_verbose_when_both_false(self):
        """When both verbose and quiet are False, default to INFO."""
        logger = setup_logging(verbose=False, quiet=False)
        assert logger.level == logging.INFO


class TestGetLogger:
    """Tests for get_logger function."""

    def test_returns_logger_instance(self):
        """get_logger should return a Logger instance."""
        logger = get_logger()
        assert isinstance(logger, logging.Logger)

    def test_returns_same_logger(self):
        """get_logger should always return the same logger instance."""
        logger1 = get_logger()
        logger2 = get_logger()
        assert logger1 is logger2


class TestLoggingHelpers:
    """Tests for logging helper functions."""

    def test_log_debug_does_not_raise(self):
        """log_debug should not raise an exception."""
        setup_logging(verbose=True)
        log_debug("Test debug message")

    def test_log_info_does_not_raise(self):
        """log_info should not raise an exception."""
        setup_logging()
        log_info("Test info message")

    def test_log_warning_does_not_raise(self):
        """log_warning should not raise an exception."""
        setup_logging()
        log_warning("Test warning message")

    def test_log_error_does_not_raise(self):
        """log_error should not raise an exception."""
        setup_logging()
        log_error("Test error message")
