#!/usr/bin/env python3
import sys
import json
from serpapi import SerpApiClient

params = {
  "api_key": "de042f985eeb4a064aae150e4dea063d8eb88309c008a32890a24560785e1ed3",
  "engine": "google_scholar",
  "q": "cretan, OR crete, OR kriti",
  "hl": "en",
  #  "num": "9"
}

search = SerpApiClient(params)
with open('results/scholar_crete.tsv', 'w') as output:
    
    output.write("page_number" + "\t" + "position" + "\t" + "result_type" + "\t" + "title"+"\t" + "link" + "\t" + "result_id" + "\t" + "publication_info_summary" + "\t" + "snippet" + "\n")

    for page in search.pagination():
        page_number = page.get('serpapi_pagination', {}).get('current')
        #print(f"Currently extracting page #{page_number}..")
    
        for result in page.get("organic_results", []):
            position = result["position"]
            title = result["title"]
            publication_info_summary = result["publication_info"]["summary"]
            result_id = result["result_id"]
            link = result.get("link")
            result_type = result.get("type")
            snippet = result.get("snippet")
    
            output.write(str(page_number) + "\t" + str(position+1) + "\t" + str(result_type) + "\t" + str(title)+"\t" + str(link) + "\t" + str(result_id) + "\t" + str(publication_info_summary) + "\t" + str(snippet) + "\n")
#            publications.append({
#                "page_number": page_number,
#                "position": position + 1,
#                "result_type": result_type,
#                "title": title,
#                "link": link,
#                "result_id": result_id,
#                "publication_info_summary": publication_info_summary,
#                "snippet": snippet,
#            })
    
    
