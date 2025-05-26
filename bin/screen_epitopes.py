#!/usr/bin/env python3

# screen_epitopes.py
import argparse
import time
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def load_epitopes_from_csv(file):
    df = pd.read_csv(file)
    return df['sequence'].dropna().astype(str).unique().tolist()


def submit_to_allertop(sequence, driver):
    driver.get("https://www.ddgpharmfac.net/allertop_test/")


    # Add fallback delay in case the site is slow
    time.sleep(3)

    # Updated to new input field
    WebDriverWait(driver, 40).until(
        EC.presence_of_element_located((By.ID, "id_protein"))
    )
    textarea = driver.find_element(By.ID, "id_protein")
    textarea.clear()
    textarea.send_keys(sequence)

    # Click the Submit button
    submit_button = driver.find_element(By.XPATH, "//button[@type='submit']")
    submit_button.click()

    # Wait for the result
    WebDriverWait(driver, 40).until(
        EC.presence_of_element_located((By.XPATH, "//pre"))
    )
    result = driver.find_element(By.XPATH, "//pre").text
    return result.strip()


def submit_to_toxinpred(epitope, driver):
    driver.get("https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.html")
    textarea = driver.find_element(By.NAME, "sequence")
    textarea.clear()
    textarea.send_keys(epitope)
    driver.find_element(By.NAME, "submit").click()

    try:
        result_table = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.XPATH, "//table[contains(.,'Peptide')]"))
        )
        rows = result_table.find_elements(By.TAG_NAME, "tr")
        for row in rows:
            if epitope in row.text:
                return "Non-Toxic" if "non-toxic" in row.text.lower() else "Toxic"
        return "Not Found"
    except:
        return "Error"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bcell")
    parser.add_argument("--tcelli")
    parser.add_argument("--tcellii")
    parser.add_argument("--output")
    args = parser.parse_args()

    if not args.output:
        print("Output file is required.")
        return

    epitopes = set()
    for file in [args.bcell, args.tcelli, args.tcellii]:
        epitopes.update(load_epitopes_from_csv(file))

    options = Options()
    options.add_argument("--headless")
    driver = webdriver.Chrome(options=options)

    results = []
    for ep in epitopes:
        allergen = submit_to_allertop(ep, driver)
        time.sleep(2)
        toxin = submit_to_toxinpred(ep, driver)
        time.sleep(2)
        results.append({"Epitope": ep, "Allergenicity": allergen, "Toxicity": toxin})

    driver.quit()
    pd.DataFrame(results).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
