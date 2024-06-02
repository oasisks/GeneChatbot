import os

from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv

"""
Simple wrapper to get the mongoDB
"""
load_dotenv()

uri = os.getenv("MONGO_URI")

mongo_client = MongoClient(uri)

try:
    mongo_client.admin.command('ping')
    print("Pinged your deployment. You successfully connected to MongoDB!")
except Exception as e:
    print(e)
