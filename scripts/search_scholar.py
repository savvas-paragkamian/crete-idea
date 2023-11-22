#!/usr/bin/env python3
import sys
import json
from serpapi import GoogleSearch

params = {
  "api_key": "de042f985eeb4a064aae150e4dea063d8eb88309c008a32890a24560785e1ed3",
  "engine": "google_scholar",
  "q": "cretan, OR crete, OR kriti",
  "hl": "en",
  "num": "914"
}

search = GoogleSearch(params)
results = search.get_dict()


with open('results/scholar_crete.json', 'w') as fp:
    json.dump(results, fp)
