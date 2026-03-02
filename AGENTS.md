# AGENTS.md

## Build/Test Commands

This repository contains multiple Docker images for bioinformatics environments. Each project has its own subdirectory with a Dockerfile.

### Available Projects
- `circosplot` - Circos visualization environment (micromamba-based)
- `crc-rna-atlas` - CRC RNA atlas analysis environment
- `grf-2025` - Conda-based GRF environment with Circos
- `moonfinder3` - R/Python bioinformatics environment
- `mpage` - MPage analysis environment
- `mpage-paper` - MPage paper environment
- `rna-modification-signature-paper` - RNA modification signature analysis
- `snakemake` - Snakemake workflow environment

### Docker Build (Single Project)
```bash
# Build any project locally (replace <project> with project name)
docker build -t <project>:latest ./<project>/

# Example: Build moonfinder3
docker build -t moonfinder3:latest ./moonfinder3/

# Build with specific platform
docker build --platform linux/amd64 -t moonfinder3:latest ./moonfinder3/

# Build with no cache
docker build --no-cache -t moonfinder3:latest ./moonfinder3/
```

### Test Docker Image
```bash
# Run container interactively
docker run -it --rm <project>:latest bash

# Test R packages (R-based images)
docker run --rm moonfinder3:latest Rscript -e "library(ggplot2); print('OK')"

# Test Python packages
docker run --rm moonfinder3:latest python -c "import pandas; print('OK')"

# Test specific package
docker run --rm moonfinder3:latest Rscript -e "library(Seurat); print(packageVersion('Seurat'))"
```

### CI/CD Workflow
GitHub Actions builds multi-arch images (amd64, arm64) and pushes to Docker Hub. Triggered on:
- Push to `<project>/**` paths
- Releases
- Manual dispatch (`workflow_dispatch`)

## Code Style Guidelines

### Dockerfile
- Use specific version tags for base images (e.g., `btrspg/vscode-base:0.0.6`)
- Order: `FROM`, `ARG`, `ENV`, `ADD`/`COPY`, then `RUN`
- Chain `RUN` commands with `&&` to reduce layers
- Clean up package manager caches in same `RUN` layer
- Add files to `/tmp/` for installation, remove with `rm -rf /tmp/*`
- Use `USER` directive for non-root execution

### R Scripts
- Use 4-space indentation
- Wrap installations in `tryCatch()` blocks
- Exit with status 1 on errors: `quit(status = 1)`
- Use `message()` for error reporting
- Load packages with `library()` or `require()`
- Use `file.path()` for cross-platform paths

### Bash Scripts
- Use `set -e` to exit on errors
- Quote variables: `"${VAR}"`
- Use functions for reusable code blocks
- Add comments for complex operations

### Package Lists (txt files)
- One package name per line
- No trailing whitespace
- Alphabetical ordering preferred
- Empty line at end of file

### GitHub Actions YAML
- Use descriptive job names
- Pin action versions (e.g., `actions/checkout@v3`)
- Use `env` for registry images
- Structure: `name`, `on`, `env`, `jobs`
- Cache Docker layers with `cache-from`/`cache-to`

### Naming Conventions
- Files: lowercase with underscores (e.g., `install_packages.R`)
- Directories: lowercase with hyphens (e.g., `grf-2025`)
- Docker images: lowercase with hyphens
- Workflow files: match project directory name
- Environment variables: `UPPERCASE_WITH_UNDERSCORES`

### Error Handling
- Docker: Chain commands with `&&` or use `set -e`
- R: Always use `tryCatch()` for installations
- Bash: Use `set -euo pipefail` for strict mode
- CI: Use `if-no-files-found: error` for artifacts

### Security
- Never commit secrets or credentials
- Use GitHub Secrets (`secrets.DH_USER`, `secrets.DH_TOKEN`)
- Pin base image versions, avoid `latest` tag
- Use non-root `USER` when possible

### Repository Structure
```
.
├── <project>/
│   ├── Dockerfile              # Main build file
│   ├── install_packages.R      # R package installer
│   ├── CRAN_packages.txt       # CRAN packages
│   ├── bioconductor_packages.txt
│   ├── devtools_packages.txt   # GitHub packages
│   └── requirements.txt        # Python packages
├── .github/workflows/          # CI/CD configs
├── daemon.json                 # Docker daemon config
└── AGENTS.md                   # This file
```

### Pre-commit Checklist
- [ ] Dockerfile builds successfully (`docker build`)
- [ ] Package lists have no duplicates
- [ ] YAML syntax valid
- [ ] All text files end with newline
- [ ] No trailing whitespace
