#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Logging configuration for diphtOscan.

This module provides centralized logging configuration with support for
verbose and quiet modes.
"""

import logging
import sys
from typing import Optional


def setup_logging(verbose: bool = False, quiet: bool = False) -> logging.Logger:
    """
    Configure logging for diphtOscan based on verbosity settings.

    Parameters
    ----------
    verbose : bool, optional
        Enable debug-level logging with timestamps. Default is False.
    quiet : bool, optional
        Suppress all output except errors. Default is False.

    Returns
    -------
    logging.Logger
        Configured logger instance for diphtOscan.

    Notes
    -----
    Logging levels:
    - quiet mode: ERROR only
    - normal mode: INFO
    - verbose mode: DEBUG with timestamps
    """
    logger = logging.getLogger('diphtoscan')
    logger.handlers.clear()
    logger.propagate = False

    handler = logging.StreamHandler(sys.stderr)

    if quiet:
        logger.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
    elif verbose:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        )
    else:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')

    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def get_logger() -> logging.Logger:
    """
    Get the diphtOscan logger instance.

    Returns
    -------
    logging.Logger
        The diphtOscan logger. If setup_logging() hasn't been called,
        returns an unconfigured logger.
    """
    return logging.getLogger('diphtoscan')


def log_debug(message: str) -> None:
    """
    Log a debug message.

    Parameters
    ----------
    message : str
        Message to log at DEBUG level.
    """
    get_logger().debug(message)


def log_info(message: str) -> None:
    """
    Log an info message.

    Parameters
    ----------
    message : str
        Message to log at INFO level.
    """
    get_logger().info(message)


def log_warning(message: str) -> None:
    """
    Log a warning message.

    Parameters
    ----------
    message : str
        Message to log at WARNING level.
    """
    get_logger().warning(message)


def log_error(message: str) -> None:
    """
    Log an error message.

    Parameters
    ----------
    message : str
        Message to log at ERROR level.
    """
    get_logger().error(message)
