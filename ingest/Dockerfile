FROM python:3.7

WORKDIR /code

# Upgrade pip
RUN pip install --upgrade pip

# copy the dependencies file to the working directory
COPY code/requirements.txt .
COPY code/app.py .

# install dependencies
RUN pip install -r requirements.txt

# command to run on container start
CMD [ "python", "app.py" ]
