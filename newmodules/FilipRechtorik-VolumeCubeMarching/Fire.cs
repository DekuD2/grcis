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

namespace FilipRechtorik
{
  [Serializable]
  public class Fire : DefaultSceneNode, ISolid
  {
    public enum Mode { Fire, Solid };

    public Mode RenderMode = Mode.Fire;

    static Stopwatch distSw = new Stopwatch();

    float[,,] volume;

    // Size in each direction
    public int Sx { get; }
    public int Sy { get; }
    public int Sz { get; }

    // Inverse size (to avoid division later on)
    public float ISx { get; } // 1 / (Sx - 1) ...
    public float ISy { get; }
    public float ISz { get; }

    int maxDim;

    // Determines the value treshold that creates the sphere
    public float Threshold = 0.5f;

    // Determines the radius to which the (value - threshold) is linearly converted into
    public float MinRadius = 0.1f;
    public float MaxRadius = 0.5f;

    double epsilon;

    public Fire (int sx, int sy, int sz)
    {
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

    public void CloudsRandomize (int numPoints = 15)
    {
      Random rnd = new Random();

      Vector3[] points = new Vector3[numPoints];
      for (int i = 0; i < numPoints; i++)
        points[i] = new Vector3((float)rnd.NextDouble(), (float)rnd.NextDouble(), (float)rnd.NextDouble());

      for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++)
          for (int k = 0; k < Sz; k++)
          {
            Vector3 curr = new Vector3(i * ISx, j * ISy, k * ISz);

            // Find the distance to the closest point
            float min = float.MaxValue;
            for (int p = 0; p < numPoints; p++)
            {
              float dist = (curr - points[p]).LengthSquared;
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

    class CellRayIntersection
    {
      public int Cx;
      public int Cy;
      public int Cz;
      public double MinT;
      public double MaxT;
    }

    IEnumerable<CellRayIntersection> AllCells ()
    {
      for (int i = 0; i < (Sx - 1); i++)
        for (int j = 0; j < (Sx - 1); j++)
          for (int k = 0; k < (Sx - 1); k++)
            yield return new CellRayIntersection()
            {
              Cx = i,
              Cy = j,
              Cz = k,
              MinT = 0,
              MaxT = 50000
            };
      yield break;
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

    /// <param name="gridPos">Ray position in the grid space, where (0,0,0) and (cellsX, cellsY, cellsZ) are corners the grid.</param>
    /// <param name="gridDir">Ray direction in the grid space.</param>
    IEnumerable<CellRayIntersection> GridRayIntersectingCells (Vector3d p0, Vector3d p1)
    {
      Vector3d gridPos = new Vector3d((p0.X + 0.5) * (Sx - 1), (p0.Y + 0.5) * (Sy - 1), (p0.Z + 0.5) * (Sz - 1));
      Vector3d gridDir = new Vector3d(p1.X * (Sx - 1), p1.Y * (Sy - 1), p1.Z * (Sz - 1));

      double dx = 0; // distance of x to get into the cube "X" region (between X = -0.5 and X = 0.5 parallel planes)
      double dy = 0;
      double dz = 0;

      //Debug.WriteLine("FROM " + gridPos + " and dir " + gridDir);

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

      gridPos += t * gridDir;

      //// Ray does not intersect the cube at all
      //if (gridPos.X > cellsX || gridPos.X < 0 || gridPos.Y > cellsY || gridPos.Y < 0 || gridPos.Z > cellsZ || gridPos.Z < 0)
      //  yield break;

      // How much t to pass one grid cell
      Vector3d tPerCell = new Vector3d(1/Math.Abs(gridDir.X), 1/Math.Abs(gridDir.Y), 1/Math.Abs(gridDir.Z));
      // T remaining to next cell in each direction
      Vector3d tToNext = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
      // Ray direction
      int dirX = Math.Sign(gridDir.X);
      int dirY = Math.Sign(gridDir.Y);
      int dirZ = Math.Sign(gridDir.Z);

      // Plane that the ray passed by last.
      int lastPlaneX = dirX >= 0 ? (int)gridPos.X : (int)Math.Ceiling(gridPos.X);
      int lastPlaneY = dirY >= 0 ? (int)gridPos.Y : (int)Math.Ceiling(gridPos.Y);
      int lastPlaneZ = dirZ >= 0 ? (int)gridPos.Z : (int)Math.Ceiling(gridPos.Z);

      // Calculate the original tToNext. Leave infinity for unattainable directions.
      if (dirX != 0)
        tToNext.X = ((lastPlaneX + dirX) - gridPos.X) / gridDir.X;
      if (dirY != 0)
        tToNext.Y = ((lastPlaneY + dirY) - gridPos.Y) / gridDir.Y;
      if (dirZ != 0)
        tToNext.Z = ((lastPlaneZ + dirZ) - gridPos.Z) / gridDir.Z;

      // Calculate the current cell position
      int cellX = Math.Min((int)(lastPlaneX + 0.5 * dirX), (Sx - 1) - 1);
      int cellY = Math.Min((int)(lastPlaneY + 0.5 * dirY), (Sy - 1) - 1);
      int cellZ = Math.Min((int)(lastPlaneZ + 0.5 * dirZ), (Sz - 1) - 1);

      bool Inside ()
        => cellX >= 0 && cellX < (Sx - 1) &&
           cellY >= 0 && cellY < (Sy - 1) &&
           cellZ >= 0 && cellZ < (Sz - 1);

      while (Inside())
      {
        double dt = Math.Min(Math.Min(tToNext.X, tToNext.Y), tToNext.Z);

        yield return new CellRayIntersection()
        {
          Cx = cellX,
          Cy = cellY,
          Cz = cellZ,
          MinT = t,
          MaxT = t + dt
        };

        tToNext.X -= dt;
        tToNext.Y -= dt;
        tToNext.Z -= dt;

        // If we crossed into next cell
        if (tToNext.X <= 0)
        {
          tToNext.X += tPerCell.X;
          cellX += dirX;
        }
        if (tToNext.Y <= 0)
        {
          tToNext.Y += tPerCell.Y;
          cellY += dirY;
        }
        if (tToNext.Z <= 0)
        {
          tToNext.Z += tPerCell.Z;
          cellZ += dirZ;
        }

        t += dt;
      }
    }

    bool IntersectCellSolid (LinkedList<Intersection> result, Vector3d p0, Vector3d p1, CellRayIntersection cell)
    {
      double SphereDist (Vector3d point, Vector3d corner, double radius)
      {
        Vector3d diff = point - corner;

        // Firstly, stretch the relative distance (whole grid), so that the grid has size 1
        diff.X *= (Sx - 1);
        diff.Y *= (Sy - 1);
        diff.Z *= (Sz - 1);

        // Take the length of that whole distance and subtract the radius of now sphere
        var len = diff.Length - radius;

        // Scale the vector by len, so we get the vector from edge of the sphere to the point
        diff *= len;

        // Convert back to spheroid
        diff.X *= ISx;
        diff.Y *= ISy;
        diff.Z *= ISz;

        // Return the actual length of the vector from edge of the sphere to the point
        return diff.Length;
      }

      double? CornerDist (Vector3d point, int cx, int cy, int cz)
      {
        // Check that the corner is valid
        if (cx < 0 || cy < 0 || cz < 0 || cx >= Sx || cy >= Sy || cz >= Sz)
          return null;

        var r = volume[cx, cy, cz];
        r -= Threshold;
        if (r > 0)
        {
          Vector3d corner = new Vector3d(-0.5) + new Vector3d(cx * ISx, cy * ISy, cz * ISz);
          double radius = MinRadius + (MaxRadius - MinRadius) * r / (1 - Threshold);
          return SphereDist(point, corner, radius);
        }
        else
          return null;
      }

      double InterpolateDistances (List<double> distances)
      {
        double k = 50;
        double res = -Math.Log(distances.Select(d => Math.Exp(-k * d)).Sum()) / k;
        return res;

        return distances.Min();
      }

      // blendRange of 1 = just the cube around the point.
      // blend range of 2 = 9 cubes around the point
      double Dist (Vector3d point, int blendRange = 2)
      {
        List<double> distances = new List<double>();
        for (int i = 1 - blendRange; i <= blendRange; i++)
          for (int j = 1 - blendRange; j <= blendRange; j++)
            for (int k = 1 - blendRange; k <= blendRange; k++)
            {
              double? dist = CornerDist(point, cell.Cx + i, cell.Cy + j, cell.Cz + k);
              if (dist != null)
                distances.Add(dist.Value);
            }

        if (distances.Count == 0)
          return 10; // 10 should skip the cell 100% of the time

        return InterpolateDistances(distances);
      }

      Vector3d ApproximateNormal (Vector3d point)
      {
        // Use quicker distance (less blending) for normals, since I have to calculate this value 6 times...
        var res = new Vector3d(Dist(new Vector3d(point.X + epsilon, point.Y, point.Z)) - Dist(new Vector3d(point.X - epsilon, point.Y, point.Z)),
                            Dist(new Vector3d(point.X, point.Y + epsilon, point.Z)) - Dist(new Vector3d(point.X, point.Y - epsilon, point.Z)),
                            Dist(new Vector3d(point.X, point.Y, point.Z + epsilon)) - Dist(new Vector3d(point.X, point.Y, point.Z - epsilon)));
        return res;
      }

      double t = cell.MinT;
      while (t < cell.MaxT)
      {
        // Calculate the distance to scene
        var currPoint = p0 + t * p1;
        var dist = Dist(currPoint);

        // If we are close enough, register it as a collision
        if (dist < epsilon)
        {
          var intersection = new Intersection(this)
          {
            T = t,
            Enter = true,
            Front = true,
            NormalLocal = ApproximateNormal(currPoint),
            CoordLocal = currPoint,
            Material = new PhongMaterial(new double[]{ 1, 0, 0 }, 0.3, 0.4, 0.3, 128), //currPoint.X + 0.5, currPoint.Y + 0.5, currPoint.Z + 0.5
            SurfaceColor = new double[]{ currPoint.X + 0.5, currPoint.Y + 0.5, currPoint.Z + 0.5 }
          };

          result.AddLast(intersection);
          return true;
        }

        // Adjust current poisition by adjusting t
        t += dist;
      }

      return false;
    }

    void IntersectFire (LinkedList<Intersection> result, Vector3d p0, Vector3d p1)
    {
      // How many lines?
      double t = FindStart(p0, p1);
      double step = 1d / (maxDim * p1.Length);
      double stepSize = step * p1.Length;

      bool IsInside (Vector3d point)
      {
        return
          point.X >= -0.5 && point.Y >= -0.5 && point.Z >= -0.5 &&
          point.X <= 0.5 && point.Y <= 0.5 && point.Z <= 0.5;
      }

      double[] flame0 = new double[]{ 89 / 255d, 22 / 255d, 3 / 255d };
      double[] flame1 = new double[]{ 250 / 255d, 236 / 255d, 165 / 255d };

      double Lerp (double from, double to, double k)
        => (1 - k) * from + k * to;

      void FlameIntersection (Vector3d point, double t2, double val)
      {
        val = 0.5;
        var intersection = new Intersection(this)
        {
          T = t,
          Enter = true,
          Front = true,
          NormalLocal = -p1,
          CoordLocal = point,
          Material = new PhongMaterial(new double[]{ 1, 0, 0 }, 0.3, 0.4, 0.1, 128){ Kt = 0.1 }, //currPoint.X + 0.5, currPoint.Y + 0.5, currPoint.Z + 0.5
          SurfaceColor = new double[]{ Lerp(flame0[0], flame1[0], val), Lerp(flame0[1], flame1[1], val), Lerp(flame0[2], flame1[2], val) }
        };

        result.AddLast(intersection);
      }

      Vector3d p = p0 + p1 * t;
      bool first = true; // Just in case of numerical errors
      while (IsInside(p))
      {
        first = false;

        double val = Sample(p.X + 0.5, p.Y + 0.5, p.Z + 0.5);
        FlameIntersection(p, t, val);

        t += step;
        p = p0 + p1 * t;
      }
    }

    public override LinkedList<Intersection> Intersect (Vector3d p0, Vector3d p1)
    {
      LinkedList<Intersection> result = new LinkedList<Intersection>();

      if (RenderMode == Mode.Solid)
      {
        foreach (var cell in GridRayIntersectingCells(p0, p1))
          if (IntersectCellSolid(result, p0, p1, cell))
            return result;
      }
      else
      {
        IntersectFire(result, p0, p1);
        return result;
      }
      //Console.WriteLine(distSw.Elapsed.TotalSeconds);
      return new LinkedList<Intersection>(result.ToList().OrderBy(x => x.T));
    }

    /// <param name="s0">sphere position</param>
    /// <param name="r">sphere radius</param>
    /// <returns></returns>
    public bool IntersectSphere (Vector3d s0, double r, Vector3d p0, Vector3d p1, out double t0, out double t1)
    {
      t0 = 0;
      t1 = 0;

      p0 = p0 - s0;

      double OD;
      Vector3d.Dot(ref p0, ref p1, out OD);
      double DD;
      Vector3d.Dot(ref p1, ref p1, out DD);
      double OO;
      Vector3d.Dot(ref p0, ref p0, out OO);
      double d = OD * OD + DD * (r*r - OO); // discriminant
      if (d <= 0.0)
        return false;            // no intersections

      d = Math.Sqrt(d);

      // there will be two intersections: (-OD - d) / DD, (-OD + d) / DD
      LinkedList<Intersection> result = new LinkedList<Intersection> ();
      t0 = (-OD - d) / DD;
      t1 = (-OD + d) / DD;

      return true;
    }

    /// <summary>
    /// Complete all relevant items in the given Intersection object.
    /// </summary>
    /// <param name="inter">Intersection instance to complete.</param>
    public override void CompleteIntersection (Intersection inter)
    {
      // Normal vector - no need to do anything here as NormalLocal is defined...

      // 2D texture coordinates.
      double r = Math.Sqrt(inter.CoordLocal.X * inter.CoordLocal.X + inter.CoordLocal.Y * inter.CoordLocal.Y);
      inter.TextureCoord.X = Geometry.IsZero(r)
        ? 0.0
        : (Math.Atan2(inter.CoordLocal.Y, inter.CoordLocal.X) / (2.0 * Math.PI) + 0.5);
      inter.TextureCoord.Y = Math.Atan2(r, inter.CoordLocal.Z) / Math.PI;
    }

    public void GetBoundingBox (out Vector3d corner1, out Vector3d corner2)
    {
      corner1 = new Vector3d(-1, -1, -1);
      corner2 = new Vector3d(1, 1, 1);
      //corner1 = new Vector3d(-0.5, -0.5, -0.5);
      //corner2 = new Vector3d(0.5, 0.5, 0.5);
    }

    public int getSerial ()
    {
      throw new NotImplementedException();
    }

    public object Clone ()
    {
      throw new NotImplementedException();
    }
  }
}
