project = "prodes"
root_doc = "README"
extensions: list[str] = []
exclude_patterns = [
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    "data",
    "docs",
    "scripts",
    "src",
    "tests",
]
html_theme = "alabaster"
