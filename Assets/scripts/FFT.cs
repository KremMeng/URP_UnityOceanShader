using System.Numerics;
using UnityEngine;

public class FFT
{
    public struct ComplexList
    {
        public Complex[] complices;
    }

    int N, which;
    int log2N;
    float pi2 = 2 * Mathf.PI;
    uint[] reversed;
    ComplexList[] T;
    ComplexList[] c = new ComplexList[2];

    public FFT(int N)
    {
        this.N = N;
        c[0].complices = new Complex[N];
        c[1].complices = new Complex[N];

        log2N = (int)(Mathf.Log(N) / Mathf.Log(2));
        //Debug.Log(log2N);

        reversed = new uint[N];
        for (int i = 0; i < N; i++)
        {
            reversed[i] = Reverse((uint)i);
            //Debug.Log(reversed[i]);
        }

        int pow2 = 1;
        T = new ComplexList[log2N];
        for (int i = 0; i < log2N; i++)
        {
            T[i].complices = new Complex[pow2];
            for (int j = 0; j < pow2; j++)
            {
                T[i].complices[j] = ComputeT(j, pow2 * 2);
            }
            pow2 *= 2;
        }

        which = 0;
    }

    //
    public uint Reverse(uint i)
    {
        uint res = 0;
        for (int j = 0; j < log2N; j++)
        {
            res = (res << 1) + (i & 1);
            i >>= 1;
        }
        return res;
    }

    public Complex ComputeT(int x, int n)
    {
        //Debug.Log(Mathf.Cos((float)pi2 * x / n));
        return new Complex(Mathf.Cos(pi2 * x / n), Mathf.Sin(pi2 * x / n));
    }

    public void ComputeFFT(Complex[] input, ref Complex[] output, int stride, int offset)
    {
        for (int i = 0; i < N; i++)
        {
            c[which].complices[i] = input[reversed[i] * stride + offset];
        }

        int loops = N >> 1;
        int size = 1 << 1;
        int sizeOver2 = 1;
        int w = 0;
        for (int i = 1; i <= log2N; i++)
        {
            //0,1
            which ^= 1;
            for (int j = 0; j < loops; j++)
            {
                for (int k = 0; k < sizeOver2; k++)
                {
                    c[which].complices[size * j + k] = c[which ^ 1].complices[size * j + k] +
                        c[which ^ 1].complices[size * j + sizeOver2 + k] * T[w].complices[k];
                }
                for (int k = sizeOver2; k < size; k++)
                {
                    c[which].complices[size * j + k] = c[which ^ 1].complices[size * j - sizeOver2 + k] -
                        c[which ^ 1].complices[size * j + k] * T[w].complices[k - sizeOver2];
                }
            }
            loops >>= 1;
            size <<= 1;
            sizeOver2 <<= 1;
            w++;
        }

        for (int i = 0; i < N; i++)
        {
            output[i * stride + offset] = c[which].complices[i];
        }
    }
}