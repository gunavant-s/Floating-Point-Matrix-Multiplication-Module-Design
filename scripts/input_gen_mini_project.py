import numpy as np
import sys, getopt
from struct import *
np.random.seed(123)

def float_to_hex(f):
  return hex(unpack('<I', pack('<f', f))[0])


def convert_header_to_dat(address, row_size, col_size):
  address = 0x00000000
  file_str = "// Row size: {}, Column size: {}\n".format(row_size, col_size)
  file_str += " @{:08X} ".format(0) + "{:04X}".format(row_size)+ "{:04X} \n".format(col_size)
  return file_str


def convert_matrix_to_dat(matrix, address, row_size, col_size, file_name, file_str):
  for i,row in enumerate(matrix):
    for j,elements in enumerate(row):
      file_str+="// {:.7f} \n".format(elements)
      file_str+=" @{:08X} ".format(address+(i*col_size)+j) + float_to_hex(elements).replace('0x','') + "\n"    

  with open(file_name,"w") as F:
    F.write(file_str)


def main(argv):
  try:
    opts, args = getopt.getopt(argv,"o:i:r:c:x:y:t:",["iDir=","oDir="])
  except getopt.GetoptError as err:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

  testFileName=''
  for opt, arg in opts:
    if "-r" in opt:
      input_row_size = int(arg)
    elif "-c" in opt:
      input_col_size = int(arg)
    elif "-t" in opt:
      print(f"testFileName received = {arg}")
      testFileName = arg	
    elif "-x" in opt:
      weight_row_size = int(arg)
    elif "-y" in opt:
      weight_col_size = int(arg)
      
  input_matrix = np.single(np.random.uniform(1,-1,size=[input_row_size,input_col_size]))
  # print(f"Input = {input_matrix}, shape= {np.shape(input_matrix)}, type= {type(input_matrix)}")

  weight_matrix = np.single(np.random.uniform(1,-1,size=[weight_row_size,weight_col_size])) 
  # print(f"Weight = {weight_matrix}, shape= {np.shape(weight_matrix)}, type= {type(weight_matrix)}")

  weight_matrix_transpose = np.array(weight_matrix).T.tolist()
  # print(f"Weight transpose = {weight_matrix_transpose}, shape= {np.shape(weight_matrix_transpose)}, type= {type(weight_matrix_transpose)}")

  result_matrix = np.matmul(input_matrix, weight_matrix)
  # print(f"Result = {result_matrix}, shape= {np.shape(result_matrix)}, type= {type(result_matrix)}")

  input_address = 0x00000000
  input_str = "// Row size: {}, Column size: {}\n".format(input_row_size, input_col_size)
  input_str += " @{:08X} ".format(0) + "{:08X}".format(input_row_size)+ "{:08X} \n".format(input_col_size)
  
  input_str = convert_header_to_dat(
    address = input_address,
    row_size=input_row_size,
    col_size=input_col_size
    )
  # print(f"Input matrix generation")
  convert_matrix_to_dat(
    matrix=input_matrix, 
    address=input_address+1, 
    row_size=input_row_size, 
    col_size=input_col_size, 
    file_name=f"{testFileName}_input.dat", 
    file_str=input_str)

  # print(f"Weight matrix generation")
  weight_address = 0x00000000
  weight_str = convert_header_to_dat(
    address =weight_address,
    row_size=weight_row_size,
    col_size=weight_col_size
    )
  convert_matrix_to_dat(
    matrix=weight_matrix_transpose, # Weight transpose is matrix is stored in memory to ease the controller
    address=weight_address+1, 
    row_size=weight_col_size, # Row and column size is flipped in Weight transpose matrix 
    col_size=weight_row_size, # Row and column size is flipped in Weight transpose matrix
    file_name=f"{testFileName}_weight.dat", 
    file_str=weight_str)

  result_str = ""
  result_address = 0x00000000

  convert_matrix_to_dat(
    matrix=result_matrix, 
    address=result_address, 
    row_size=input_row_size, 
    col_size=weight_col_size, 
    file_name=f"{testFileName}_result.dat", 
    file_str=result_str)


if __name__ == "__main__":
     print(sys.argv[1:])
     main(sys.argv[1:])
