"""Process-wide re-entrant lock protecting non-thread-safe OpenBabel C++ calls under Gunicorn gthread."""
import threading

OPENBABEL_LOCK = threading.RLock()
