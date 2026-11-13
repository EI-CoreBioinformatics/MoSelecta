# MoSelecta container image — installs a pinned release from PyPI.
#
ARG PYTHON_VERSION=3.13
FROM python:${PYTHON_VERSION}-slim

ARG MOSELECTA_VERSION=0.3.0
ARG PIP_INDEX_URL=https://pypi.org/simple/

LABEL org.opencontainers.image.title="MoSelecta" \
      org.opencontainers.image.description="CDS-structure-aware gene-model selection with optional genomic sequence, protein, and transcript evidence integration" \
      org.opencontainers.image.version="${MOSELECTA_VERSION}" \
      org.opencontainers.image.source="https://github.com/EI-CoreBioinformatics/moselecta" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN pip install --no-cache-dir --disable-pip-version-check \
        --index-url "${PIP_INDEX_URL}" \
        "moselecta==${MOSELECTA_VERSION}" \
    && moselecta --version \
    && moselecta-compare --help > /dev/null

WORKDIR /data
CMD ["moselecta", "--help"]
