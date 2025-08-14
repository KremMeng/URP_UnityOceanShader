using System;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using Random = System.Random;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class Ocean
{
    public struct VertexOcean
    {
        public Vector3[] vertices;// vertex
        public Vector3[] normals;// normal
        public Vector2[] uvs;//UV
        public Vector2[] htildes;// htilde0
        public Vector2[] conjugates;// htilde0mk conjugate
        public Vector3[] originals;// original position
        public Vector2[] uvs2;
    }

    public struct ComplexVectorNormal
    {
        // structure used with discrete fourier transform
        public Complex h;// wave height
        public Vector2 D;// displacement
        public Vector3 n;// normal
    }

    float g = 9.81f;                                // gravity constant
    int N, Nplus1;                          // dimension -- N should be a power of 2
    float A;                                // phillips spectrum parameter -- affects heights of waves
    Vector2 w;                              // wind parameter
    float length;                           // length parameter
    Random rand = new Random();

    Complex[] tilde, tildeSlopeX, tildeSlopeZ, tildeDx, tildeDz;
    FFT fft;                              // fast fourier transform

    public VertexOcean oceanData;                    // vertices for vertex buffer object

    public List<int> indices;                  // indicies for vertex buffer object

    int debug = 0;


    //M=N,L=Lx=Lz
    public Ocean(int N, float A, Vector2 w, float length)
    {
        tilde = new Complex[N * N];
        tildeSlopeX = new Complex[N * N];
        tildeSlopeZ = new Complex[N * N];
        tildeDx = new Complex[N * N];
        tildeDz = new Complex[N * N];

        fft = new FFT(N);

        indices = new List<int>();

        this.N = N;
        Nplus1 = N + 1;
        this.A = A;
        this.w = w;
        this.length = length;

        oceanData = new VertexOcean();
        oceanData.vertices = new Vector3[Nplus1 * Nplus1];
        oceanData.normals = new Vector3[Nplus1 * Nplus1];
        oceanData.uvs = new Vector2[Nplus1 * Nplus1];
        oceanData.htildes = new Vector2[Nplus1 * Nplus1];
        oceanData.conjugates = new Vector2[Nplus1 * Nplus1];
        oceanData.originals = new Vector3[Nplus1 * Nplus1];
        oceanData.uvs2 = new Vector2[Nplus1 * Nplus1];

        int index;
        Complex htilde0, htilde0mk_conj;
        for (int m = 0; m < Nplus1; m++)
        {
            for (int n = 0; n < Nplus1; n++)
            {
                index = m * Nplus1 + n;

                htilde0 = hTilde_0(n, m);
                htilde0mk_conj = Complex.Conjugate(hTilde_0(-n, -m));

                oceanData.htildes[index].x = (float)htilde0.Real;
                oceanData.htildes[index].y = (float)htilde0.Imaginary;
                oceanData.conjugates[index].x = (float)htilde0mk_conj.Real;
                oceanData.conjugates[index].y = (float)htilde0mk_conj.Imaginary;

                oceanData.originals[index].x = oceanData.vertices[index].x = (n - N / 2.0f) * length / N;
                oceanData.originals[index].y = oceanData.vertices[index].y = 0;
                oceanData.originals[index].z = oceanData.vertices[index].z = (m - N / 2.0f) * length / N;

                oceanData.normals[index].x = 0.0f;
                oceanData.normals[index].y = 1.0f;
                oceanData.normals[index].z = 0.0f;

                oceanData.uvs[index].x = (float)n / N;
                oceanData.uvs[index].y = (float)m / N;
            }
        }

        for (int m = 0; m < N; m++)
        {
            for (int n = 0; n < N; n++)
            {
                index = m * Nplus1 + n;

                indices.Add(index);
                indices.Add(index + Nplus1);
                indices.Add(index + Nplus1 + 1);
                indices.Add(index);
                indices.Add(index + Nplus1 + 1);
                indices.Add(index + 1);
            }
        }
    }

    //????????
    public Complex GaussRandomVariable()
    {
        float s = 0, u = 0, v = 0;
        while (s > 1 || s == 0)
        {
            u = (float)rand.NextDouble() * 2 - 1;
            v = (float)rand.NextDouble() * 2 - 1;

            s = u * u + v * v;
        }

        s = Mathf.Sqrt(-2 * Mathf.Log(s) / s);
        return new Complex(u * s, v * s);
    }

    //?????????
    public float Phillips(int n, int m)
    {
        Vector2 k = new Vector2(Mathf.PI * (2 * n - N) / length, Mathf.PI * (2 * m - N) / length);
        float k_Length = k.magnitude;

        if (k_Length < 0.000001) return 0.0f;

        float k_Length2 = k_Length * k_Length;
        float k_Length4 = k_Length2 * k_Length2;

        float KdotW = Vector2.Dot(k.normalized, w.normalized);
        float KdotW2 = KdotW * KdotW;

        float w_Length = w.magnitude;
        float L = w_Length * w_Length / g;
        float L2 = L * L;

        float damping = 0.001f;
        float l2 = L2 * damping * damping;

        //???????????? Mathf.Exp(-k_Length * l2)?
        float phillips = A * Mathf.Exp(-1.0f / (k_Length2 * L2)) / k_Length4 * KdotW2
                * Mathf.Exp(-k_Length2 * l2);
        return phillips;
    }

    public Complex hTilde_0(int n, int m)
    {
        Complex r = GaussRandomVariable();
        return r * Mathf.Sqrt(Phillips(n, m) / 2.0f);
    }

    Complex hTilde(float t, int n, int m)
    {
        int index = m * Nplus1 + n;

        Complex htilde0 = new Complex(oceanData.htildes[index].x, oceanData.htildes[index].y);
        Complex htilde0mk_conj = new Complex(oceanData.conjugates[index].x, oceanData.conjugates[index].y);

        float omegat = Dispersion(n, m) * t;

        float cos = Mathf.Cos(omegat);
        float sin = Mathf.Sin(omegat);

        Complex c0 = new Complex(cos, sin);
        Complex c1 = new Complex(cos, -sin);

        Complex res = htilde0 * c0 + htilde0mk_conj * c1;

        return res;
    }

    public float Dispersion(int n, int m)
    {
        float w0 = 2.0f * Mathf.PI / 200.0f;
        float kx = Mathf.PI * (2 * n - N) / length;
        float kz = Mathf.PI * (2 * m - N) / length;
        
        return Mathf.Floor(Mathf.Sqrt(g * Mathf.Sqrt(kx * kx + kz * kz) / w0)) * w0;
    }

    public ComplexVectorNormal ComputerHDN(Vector2 x, float t)
    {
        Complex h = new Complex();
        Vector2 D = Vector2.zero;
        Vector3 normal = Vector3.zero;

        Complex c, htilde;
        Vector2 k;
        float kx, kz, k_Length, KdotX;

        for (int m = 0; m < N; m++)
        {
            kz = 2.0f * Mathf.PI * (m - N / 2.0f) / length;
            for (int n = 0; n < N; n++)
            {
                kx = 2.0f * Mathf.PI * (n - N / 2.0f) / length;
                k = new Vector2(kx, kz);

                k_Length = k.magnitude;
                KdotX = Vector2.Dot(k, x);

                c = new Complex(Mathf.Cos(KdotX), Mathf.Sin(KdotX));
                htilde = hTilde(t, n, m) * c;

                h += htilde;

                normal += new Vector3(-kx * (float)htilde.Imaginary, 0.0f,
                    -kz * (float)htilde.Imaginary);

                if (k_Length < 0.000001) continue;
                D += new Vector2(kx / k_Length * (float)htilde.Imaginary, kz / k_Length *
                    (float)htilde.Imaginary);
            }
        }

        normal = new Vector3(0f, 1f, 0f) - normal;
        ComplexVectorNormal complex = new ComplexVectorNormal();
        complex.h = h;
        complex.D = D;
        complex.n = normal.normalized;
        return complex;
    }

    public void EvaluateWaves(float t)
    {
        float kx, kz, len, lambda = -1.0f;
        int index, index1;

        for (int m = 0; m < N; m++)
        {
            kz = Mathf.PI * (2.0f * m - N) / length;
            for (int n = 0; n < N; n++)
            {
                kx = Mathf.PI * (2.0f * n - N) / length;
                len = Mathf.Sqrt(kx * kx + kz * kz);
                index = m * N + n;

                tilde[index] = hTilde(t, n, m);
                tildeSlopeX[index] = tilde[index] * new Complex(0, kx);
                tildeSlopeZ[index] = tilde[index] * new Complex(0, kz);

                if (len < 0.000001f)
                {
                    tildeDx[index] = new Complex(0f, 0f);
                    tildeDz[index] = new Complex(0f, 0f);
                }
                else
                {
                    tildeDx[index] = tilde[index] * new Complex(0f, -kx / len);
                    tildeDz[index] = tilde[index] * new Complex(0f, -kz / len);
                }
            }
        }

        for (int m = 0; m < N; m++)
        {
            fft.ComputeFFT(tilde, ref tilde, 1, m * N);
            fft.ComputeFFT(tildeSlopeX, ref tildeSlopeX, 1, m * N);
            fft.ComputeFFT(tildeSlopeZ, ref tildeSlopeZ, 1, m * N);
            fft.ComputeFFT(tildeDx, ref tildeDx, 1, m * N);
            fft.ComputeFFT(tildeDz, ref tildeDz, 1, m * N);
        }

        for (int n = 0; n < N; n++)
        {
            fft.ComputeFFT(tilde, ref tilde, N, n);
            fft.ComputeFFT(tildeSlopeX, ref tildeSlopeX, N, n);
            fft.ComputeFFT(tildeSlopeZ, ref tildeSlopeZ, N, n);
            fft.ComputeFFT(tildeDx, ref tildeDx, N, n);
            fft.ComputeFFT(tildeDz, ref tildeDz, N, n);
        }

        int sign;
        int[] signs = { 1, -1 };
        Vector3 normal;
        Vector3[] offsets = ComputeOffset();
        debug = 0;
        for (int m = 0; m < N; m++)
        {
            for (int n = 0; n < N; n++)
            {
                index = m * N + n;
                index1 = m * Nplus1 + n;

                //???n+m?????????sign=signs[1]=-1?????n+m????????sign=signs[0]=1
                sign = signs[(n + m) & 1];

                tilde[index] = tilde[index] * sign;

                tildeDx[index] = tildeDx[index] * sign;
                tildeDz[index] = tildeDz[index] * sign;


                Vector3 offset = new Vector3(lambda * (float)tildeDx[index].Real, (float)tilde[index].Real,
                    lambda * (float)tildeDz[index].Real);

                oceanData.vertices[index1].y = offset.y;
                oceanData.vertices[index1].x = oceanData.originals[index1].x + offset.x;
                oceanData.vertices[index1].z = oceanData.originals[index1].z + offset.z;

                // normal
                tildeSlopeX[index] = tildeSlopeX[index] * sign;
                tildeSlopeZ[index] = tildeSlopeZ[index] * sign;
                normal = new Vector3(-(float)tildeSlopeX[index].Real, 1f, -(float)tildeSlopeZ[index].Real);
                oceanData.normals[index1] = normal;

                Vector2 jaco = new Vector2(ComputeJacobian(n, m, offsets), 0);
                //DebugLog(jaco.x.ToString());
                oceanData.uvs2[index1] = jaco;

                int number;
                if (n == 0 && m == 0)
                {
                    number = index1 + N + Nplus1 * N;
                    oceanData.vertices[number].y = offset.y;
                    oceanData.vertices[number].x = oceanData.originals[number].x + offset.x;
                    oceanData.vertices[number].z = oceanData.originals[number].z + offset.z;
                    oceanData.normals[number] = normal;
                    oceanData.uvs2[number] = jaco;
                }

                if (n == 0)
                {
                    number = index1 + N;
                    oceanData.vertices[number].y = offset.y;
                    oceanData.vertices[number].x = oceanData.originals[number].x + offset.x;
                    oceanData.vertices[number].z = oceanData.originals[number].z + offset.z;
                    oceanData.normals[number] = normal;
                    oceanData.uvs2[number] = jaco;
                }

                if (m == 0)
                {
                    number = index1 + Nplus1 * N;
                    oceanData.vertices[number].y = offset.y;
                    oceanData.vertices[number].x = oceanData.originals[number].x + offset.x;
                    oceanData.vertices[number].z = oceanData.originals[number].z + offset.z;
                    oceanData.normals[number] = normal;
                    oceanData.uvs2[number] = jaco;
                }
            }
        }
    }

    void DebugLog(string a)
    {
        if (debug > 500) return;
        debug++;
        Debug.Log(a);
    }

    Vector3[] ComputeOffset()
    {
        Vector3[] offsets = new Vector3[N * N];
        float lambda = -1.0f;
        int index;
        int[] signs = { 1, -1 };
        for (int m = 0; m < N; m++)
        {
            for (int n = 0; n < N; n++)
            {
                index = m * N + n;
                int sign = signs[(n + m) & 1];
                offsets[index] = new Vector3(lambda * (float)(tildeDx[index] * sign).Real,
                    (float)(tilde[index] * sign).Real,
                    lambda * (float)(tildeDz[index] * sign).Real);
            }
        }
        return offsets;
    }

    float ComputeJacobian(int n, int m, Vector3[] offsets)
    {
        int up = m - 1 < 0 ? (N - 1) * N + n : (m - 1) * N + n;
        int down = m + 1 > N - 1 ? n : (m + 1) * N + n;
        int left = n - 1 < 0 ? m * N + (N - 1) : m * N + n - 1;
        int right = n + 1 > N - 1 ? m * N : m * N + n + 1;

        Vector3 ddx = offsets[down] - offsets[up];
        Vector3 ddz = offsets[right] - offsets[left];

        float jacobian = (1.0f + ddx.x) * (1.0f + ddz.z) - ddx.z * ddz.x;
        return jacobian;
    }
}