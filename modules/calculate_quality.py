def calculate_quality(quality: str) -> float:
"""
This function takes read quality string as input and 
returns quality score based on Phred33 algorithm
Characters are converted in ASCII with ord() builtin 
function and then 33 substracted.
See Phred33 algorithm for reference.
Function returns mean quality for string
"""
    if not quality: 
        return 0.0
    qualities = []
    for i in quality:
        qualities.append(ord(i)-33)
    return (sum(qualities))/len(qualities)