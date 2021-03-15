import boto3
import os
import tarfile

if __name__ == "__main__":
    s3 = boto3.client('s3', aws_access_key_id="AKIAY6UR252SQUQ3OSWZ",
                      aws_secret_access_key="08LQj"
                                            "+ryk9SMojG18vERXKKzhNSYk5pLhAjrIAVX")
    output_path = "./data.tar.gz"
    with open(output_path, 'wb') as f:
        s3.download_fileobj('definedproteins', "data.tar.gz", f)
    assert os.path.isfile(output_path)
    print("Download succeeded")
    tar = tarfile.open(output_path, "r:gz")
    tar.extractall()
    tar.close()
    os.remove(output_path)