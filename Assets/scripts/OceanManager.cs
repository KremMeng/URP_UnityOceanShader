using UnityEngine;

public class OceanManager : MonoBehaviour
{
    public int N = 32;
    public float WaveHeight = 1f;
    public Vector2 Wind = new Vector2(1f, 1f);
    public float Length = 1f;

    Ocean _ocean;
    Mesh _mesh;

    void Start()
    {
        _ocean = new Ocean(N, WaveHeight, Wind, Length);
        _mesh = GetComponent<MeshFilter>().mesh;
        //_mesh.Clear();
    }

    private void FixedUpdate()
    {
        //Debug.Log(Time.time);
        _ocean.EvaluateWaves(Time.time);
        _mesh.vertices = _ocean.oceanData.vertices;
        _mesh.triangles = _ocean.indices.ToArray();
        _mesh.uv = _ocean.oceanData.uvs;
        _mesh.normals = _ocean.oceanData.normals;
        _mesh.uv2 = _ocean.oceanData.uvs2;
        //_mesh.RecalculateNormals();
    }
}