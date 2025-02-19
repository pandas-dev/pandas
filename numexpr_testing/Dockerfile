# Use Python 3.9 as base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Copy requirements first to leverage Docker cache
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy test files
COPY test_evaluate_numexpr.py .

# Command to run tests
CMD ["pytest", "test_evaluate_numexpr.py", "-v"]