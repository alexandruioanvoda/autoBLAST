def reverse(a_string):
    return a_string[::-1]

def remove_dashes(a_string):
    return a_string.replace("-", "")

def complement(a_string):
    a_string = a_string.replace("A", "1")
    a_string = a_string.replace("a", "2")
    a_string = a_string.replace("C", "3")
    a_string = a_string.replace("c", "4")
    a_string = a_string.replace("G", "5")
    a_string = a_string.replace("g", "6")
    a_string = a_string.replace("T", "7")
    a_string = a_string.replace("t", "8")
    
    a_string = a_string.replace("1", "T")
    a_string = a_string.replace("2", "t")
    a_string = a_string.replace("3", "G")
    a_string = a_string.replace("4", "g")
    a_string = a_string.replace("5", "C")
    a_string = a_string.replace("6", "c")
    a_string = a_string.replace("7", "A")
    a_string = a_string.replace("8", "a")
    
    return a_string
