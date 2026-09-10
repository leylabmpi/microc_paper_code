#!/usr/bin/env bash
set -euo pipefail

RECORD_ID="22687354"
FILE_NAME="data.zip"

echo "Downloading ${FILE_NAME}..."
curl -L -o "${FILE_NAME}" "https://zenodo.org/records/${RECORD_ID}/files/${FILE_NAME}?download=1"

echo "Unzipping..."
unzip -o "${FILE_NAME}"
rm "${FILE_NAME}"

echo "Done. Folder structure restored under data/"