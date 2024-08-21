from functools import wraps
import logging
from time import perf_counter
from typing import Callable, Any, Optional


def time_it(message: Optional[str] = None) -> Callable[..., Any]:
    """Decorator to measure the execution time of a function with an optional formatted message."""

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            start_time = perf_counter()
            result = func(*args, **kwargs)
            end_time = perf_counter()
            elapsed_time = end_time - start_time

            # Combine positional and keyword arguments into a single dictionary
            arg_names = func.__code__.co_varnames[:func.__code__.co_argcount]
            arg_dict = dict(zip(arg_names, args))
            all_kwargs = {**arg_dict, **kwargs}

            func_name = func.__name__

            # If a message is provided, format it with combined arguments
            if message:
                try:
                    final_message = message.format(**all_kwargs)
                except KeyError as e:
                    raise KeyError(f"Missing key in message format: {e}")
            else:
                final_message = func.__name__

            log_message = f"Execution time for {final_message}:\t{elapsed_time:.2f} seconds\t{func_name}"
            logging.info(log_message)
            print(log_message)
            return result

        return wrapper

    return decorator
