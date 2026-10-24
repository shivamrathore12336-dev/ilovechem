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

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Run DB migration and collect static
RUN python manage.py migrate

EXPOSE 8000

CMD gunicorn chemsite.wsgi:application --bind 0.0.0.0:${PORT:-8000} --workers 3 --timeout 120
