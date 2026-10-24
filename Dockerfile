FROM python:3.11-slim

# System dependencies for OpenBabel and RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libopenbabel-dev \
    openbabel \
    libxrender1 \
    libxext6 \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set up non-root user (UID 1000) for security and Hugging Face Spaces compatibility
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH

WORKDIR $HOME/app

COPY --chown=user:user requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

COPY --chown=user:user . .

# Run DB migration for SQLite session/cache
RUN python manage.py migrate

EXPOSE 10000

CMD ["sh", "-c", "gunicorn chemsite.wsgi:application --bind 0.0.0.0:${PORT:-10000} --workers ${WEB_CONCURRENCY:-1} --timeout 120"]
