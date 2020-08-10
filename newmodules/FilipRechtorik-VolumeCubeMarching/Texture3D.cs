using MathSupport;
using OpenTK;
using Scene3D;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using Rendering;
using System.Runtime.CompilerServices;
using System.Collections.Specialized;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Threading;
using System.Windows.Forms.VisualStyles;
using System.Security.AccessControl;
using System.CodeDom.Compiler;
using OpenTK.Graphics.ES11;
using System.Security.Cryptography;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Drawing;
using System.Configuration;
using Utilities;
using System.Drawing.Drawing2D;

namespace FilipRechtorik
{
  [Serializable]
  public class Texture3D : DefaultSceneNode, ISolid
  {
    protected int seed;
    protected Random random;

    protected float[,,] volume;

    // Size in each direction
    public int Sx { get; }
    public int Sy { get; }
    public int Sz { get; }

    // Inverse size (to avoid division later on)
    public float ISx { get; } // 1 / (Sx - 1) ...
    public float ISy { get; }
    public float ISz { get; }

    protected int maxDim;
    protected double epsilon;

    public void SetSeed(int seed)
    {
      this.seed = seed;
      this.random = new Random(seed);
    }

    public Texture3D (int sx, int sy, int sz)
    {
      Random rnd = new Random();
      seed = rnd.Next();
      random = new Random(seed);

      volume = new float[sx, sy, sz];
      maxDim = Math.Max(sx, Math.Max(sy, sz));
      epsilon = 0.01 / maxDim;

      Sx = sx;
      Sy = sy;
      Sz = sz;
      ISx = 1 / (Sx - 1f);
      ISy = 1 / (Sy - 1f);
      ISz = 1 / (Sz - 1f);
    }

    // Input: coordinates 0-1
    public double Sample (double x, double y, double z)
    {
      // Scale to grid
      x *= (Sx - 1);
      y *= (Sy - 1);
      z *= (Sz - 1);

      // Lower coordinates in the texture
      int lx = (int)x;
      int ly = (int)y;
      int lz = (int)z;

      // Uper coordinates
      int ux = lx + 1;
      int uy = ly + 1;
      int uz = lz + 1;

      // Proportions from current to next value
      double px = x - lx;
      double py = y - ly;
      double pz = z - lz;

      // Inverse proportions
      double ipx = 1 - px;
      double ipy = 1 - py;
      double ipz = 1 - pz;

      double value =
        volume[lx, ly, lz] * ipx * ipy * ipz +
        volume[lx, ly, uz] * ipx * ipy * pz +
        volume[lx, uy, lz] * ipx * py * ipz +
        volume[lx, uy, uz] * ipx * py * pz +
        volume[ux, ly, lz] * px * ipy * ipz +
        volume[ux, ly, uz] * px * ipy * pz +
        volume[ux, uy, lz] * px * py * ipz +
        volume[ux, uy, uz] * px * py * pz;

      return value;
    }

    public void PointsRandomize (double min, double max, int nPoints)
    {
      for (int i = 0; i < nPoints; i++)
      {
        int x = random.Next(1, Sx - 1);
        int y = random.Next(1, Sy - 1);
        int z = random.Next(1, Sz - 1);
        double val = random.NextDouble() * (max - min) + min;
        volume[x, y, z] = (float)val;
      }
    }

    public void CloudsRandomize (int numPoints = 15)
    {
      Vector3[] points = new Vector3[numPoints];
      for (int i = 0; i < numPoints; i++)
        //points[i] = new Vector3((float)random.NextDouble(), (float)random.NextDouble(), (float)random.NextDouble());
        points[i] = new Vector3(random.Next(0, Sx), random.Next(0, Sy), random.Next(0, Sz));

      for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++)
          for (int k = 0; k < Sz; k++)
          {
            //Vector3 curr = new Vector3(i * ISx, j * ISy, k * ISz);
            Vector3 curr = new Vector3(i, j, k);

            // Find the distance to the closest point
            float min = float.MaxValue;
            for (int p = 0; p < numPoints; p++)
            {
              float dist = (curr - points[p]).LengthSquared * 0.1f / (float)maxDim;
              if (dist < min)
                min = dist;
            }

            volume[i, j, k] = 1 - (float)Math.Sqrt(Math.Sqrt(min));
          }

      ClearEdges();
    }

    public void ClearEdges ()
    {
      for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++)
          for (int k = 0; k < Sz; k++)
          {
            if (i == 0 || j == 0 || k == 0 || i == (Sx - 1) || j == (Sy - 1) || k == (Sz - 1))
              volume[i, j, k] = 0;
          }
    }

    double FindStart (Vector3d p0, Vector3d p1)
    {
      Vector3d gridPos = new Vector3d((p0.X + 0.5) * (Sx - 1), (p0.Y + 0.5) * (Sy - 1), (p0.Z + 0.5) * (Sz - 1));
      Vector3d gridDir = new Vector3d(p1.X * (Sx - 1), p1.Y * (Sy - 1), p1.Z * (Sz - 1));

      double dx = 0; // distance of x to get into the cube "X" region (between X = -0.5 and X = 0.5 parallel planes)
      double dy = 0;
      double dz = 0;

      // Add some epsilon to get inside the grid and not stay on the edge
      // Before I did this, I would end up on the edge of the cube, sometimes on the wrong side due to numerical imprecision, and that would create a noise effect.
      double e = 0.1;
      if (gridPos.X > (Sx - 1) && gridDir.X < 0)
        dx = (Sx - 1) - gridPos.X - e;
      else if (gridPos.X < 0 && gridDir.X > 0)
        dx = -gridPos.X + e;

      if (gridPos.Y > (Sy - 1) && gridDir.Y < 0)
        dy = (Sy - 1) - gridPos.Y - e;
      else if (gridPos.Y < 0 && gridDir.Y > 0)
        dy = -gridPos.Y + e;

      if (gridPos.Z > (Sz - 1) && gridDir.Z < 0)
        dz = (Sz - 1) - gridPos.Z - e;
      else if (gridPos.Z < 0 && gridDir.Z > 0)
        dz = -gridPos.Z + e;

      double t = dx / gridDir.X; // > 0 since both have same sign
      t = Math.Max(t, dy / gridDir.Y);
      t = Math.Max(t, dz / gridDir.Z);

      return t;
    }

    public void GetBoundingBox (out Vector3d corner1, out Vector3d corner2)
    {
      //corner1 = new Vector3d(-1, -1, -1);
      //corner2 = new Vector3d(1, 1, 1);
      corner1 = new Vector3d(-0.5, -0.5, -0.5);
      corner2 = new Vector3d(0.5, 0.5, 0.5);
    }

    public int getSerial ()
    {
      throw new NotImplementedException();
    }
  }
}
