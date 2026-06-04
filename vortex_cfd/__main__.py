import logging
from vortex_cfd.cli import main

# Console handler — set up once at entry point so all modules share it.
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[logging.StreamHandler()],
)

if __name__ == "__main__":
    main()
