from flask import Flask, request
from time import sleep

app = Flask(__name__)


@app.route("/run", methods=["POST"])
def hello():
  return "Hello World!"
