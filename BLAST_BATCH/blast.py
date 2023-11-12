import pandas as pd
from argparse import ArgumentParser
from tqdm import tqdm
import pickle
import torch
from transformers import GPT2LMHeadModel, GPT2Tokenizer, AutoTokenizer
import random
from selenium import webdriver
from selenium.webdriver.common.by import By
import os
import time
from tqdm import tqdm
from selenium.webdriver.firefox.options import Options
from Bio import SeqIO
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def blast(seqs):
    # os.environ['MOZ_HEADLESS'] = '1' --> works only on macOS    
    # options = Options()
    # options.add_argument('--headless')
    # driver = webdriver.Firefox(options=options)
    driver = webdriver.Firefox()
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'
    driver.get(url)
    item = driver.find_element('id','seq')
    item.send_keys(seqs)
    time.sleep(5)
    item = driver.find_element(By.CLASS_NAME,'blastbutton')
    item.click()
    item = WebDriverWait(driver, 10000, 30).until(EC.presence_of_element_located((By.ID, 'ulDnldAl'))).click()
    time.sleep(5)
    driver.find_element("link text", "Text").click()
    time.sleep(480)
    print('Finished Blast and Download')
    driver.quit()
    
    
    
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))
    