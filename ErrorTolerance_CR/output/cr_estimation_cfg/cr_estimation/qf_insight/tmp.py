import bitstring

data = [1.000000,5.000000,10.00000,20.00000,15.00000,
        1.000000,20.00000,50.00000,20.00000,5.000000,
        1.000000,15.00000,10.00000,15.00000,5.000000,
        10.00000,15.00000,10.00000,15.00000,1.000000]
f = open("data_b.dat","w")

for i in range(0,len(data)):
    num = data[i]
    # print(type(bitstring.BitArray(float=num,length=32).bin))
    f.write(bitstring.BitArray(float=num,length=32).bin)
