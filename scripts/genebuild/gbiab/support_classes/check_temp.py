import tempfile

print("Creating one named temporary file...")

temp = tempfile.NamedTemporaryFile()

try:
    print("Created file is:", temp)
    print("Name of the file is:", temp.name)
finally:
    print("Closing the temp file")
    temp.close()
