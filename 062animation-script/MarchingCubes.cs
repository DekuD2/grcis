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

namespace FilipRechtorik
{
  [Serializable]
  public class MarchingCubes : DefaultSceneNode, ISolid
  {
    //int[,,]
    // TODO: Make two bigger in each direction, let the surface be 0 everywhere
    double[,,] grid;

    // All faces (for all grid points)
    int allFacesX;
    int allFacesY;
    int allFacesZ;

    // All cells for rendering
    int allCellsX => allFacesX - 1;
    int allCellsY => allFacesY - 1;
    int allCellsZ => allFacesZ - 1;

    // Only internal faces (excluding the surface points with implicit value 0)
    public int FacesX => allFacesX - 2;
    public int FacesY => allFacesY - 2;
    public int FacesZ => allFacesZ - 2;

    // Only internal cells
    public int CellsX => allCellsX - 2;
    public int CellsY => allCellsY - 2;
    public int CellsZ => allCellsZ - 2;

    public double this[int x, int y, int z]
    {
      get => grid[x + 1, y + 1, z + 1];
      set => grid[x + 1, y + 1, z + 1] = value;
    }

    public double Start { get; set; }
    public double End { get; set; }
    public double Time { get; set; }

    public MarchingCubes (int sizeX, int sizeY, int sizeZ)
    {
      sizeX += 2;
      sizeY += 2;
      sizeZ += 2;

      this.allFacesX = sizeX;
      this.allFacesY = sizeY;
      this.allFacesZ = sizeZ;

      // Initialize randomized 3d grid
      grid = new double[sizeX, sizeY, sizeZ];
      //Random rnd = new Random();
      //for (int i = 1; i < sizeX - 1; i++)
      //  for (int j = 1; j < sizeY - 1; j++)
      //    for (int k = 1; k < sizeZ - 1; k++)
      //      grid[i, j, k] = 1.0 - Math.Sqrt(rnd.NextDouble());
      //// DEBUG: Setting specific shapes
      //grid[2, 2, 2] = 1;
      //grid[1, 2, 2] = 0.9;
      //grid[2, 1, 2] = 0.7;
      ////grid[2, 2, 1] = 0.6;
      ////grid[1, 1, 2] = 0.8;
      //grid[1, 2, 1] = 0.85;
      //grid[1, 1, 1] = 0.75;
    }

    public void Randomize ()
    {
      Random rnd = new Random();
      for (int i = 0; i < FacesX; i++)
        for (int j = 0; j < FacesY; j++)
          for (int k = 0; k < FacesZ; k++)
            this[i, j, k] = 1.0 - Math.Sqrt(rnd.NextDouble());
    }

    public void Sphere ()
    {
      Vector3d center = new Vector3d(FacesX / 2.0, FacesY / 2.0, FacesZ / 2.0);
      double radius = Math.Min(Math.Min(FacesX, FacesY), FacesZ) / 2.0;
      for (int i = 0; i < FacesX; i++)
        for (int j = 0; j < FacesY; j++)
          for (int k = 0; k < FacesZ; k++)
          {
            Vector3d curr = new Vector3d(i,j,k) - center;
            double dist = curr.Length / radius;
            if (dist < 1)
              this[i, j, k] = 1;
            else
              this[i, j, k] = 0; // 1 - dist = how far is the sample point.
          }
    }

    public void Fill ()
    {
      for (int i = 0; i < FacesX; i++)
        for (int j = 0; j < FacesY; j++)
          for (int k = 0; k < FacesZ; k++)
            this[i, j, k] = 1;
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
      for (int i = 0; i < allCellsX; i++)
        for (int j = 0; j < allCellsX; j++)
          for (int k = 0; k < allCellsX; k++)
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

    /// <param name="gridPos">Ray position in the grid space, where (0,0,0) and (cellsX, cellsY, cellsZ) are corners the grid.</param>
    /// <param name="gridDir">Ray direction in the grid space.</param>
    IEnumerable<CellRayIntersection> GridRayIntersectingCells (Vector3d gridPos, Vector3d gridDir)
    {
      double dx = 0; // distance of x to get into the cube "X" region (between X = -0.5 and X = 0.5 parallel planes)
      double dy = 0;
      double dz = 0;

      //Debug.WriteLine("FROM " + gridPos + " and dir " + gridDir);

      // Add some epsilon to get inside the grid and not stay on the edge
      // Before I did this, I would end up on the edge of the cube, sometimes on the wrong side due to numerical imprecision, and that would create a noise effect.
      double e = 0.1;
      if (gridPos.X > allCellsX && gridDir.X < 0)
        dx = allCellsX - gridPos.X - e;
      else if (gridPos.X < 0 && gridDir.X > 0)
        dx = -gridPos.X + e;

      if (gridPos.Y > allCellsY && gridDir.Y < 0)
        dy = allCellsY - gridPos.Y - e;
      else if (gridPos.Y < 0 && gridDir.Y > 0)
        dy = -gridPos.Y + e;

      if (gridPos.Z > allCellsZ && gridDir.Z < 0)
        dz = allCellsZ - gridPos.Z - e;
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
      int cellX = Math.Min((int)(lastPlaneX + 0.5 * dirX), allCellsX - 1);
      int cellY = Math.Min((int)(lastPlaneY + 0.5 * dirY), allCellsY - 1);
      int cellZ = Math.Min((int)(lastPlaneZ + 0.5 * dirZ), allCellsZ - 1);

      bool Inside ()
        => cellX >= 0 && cellX < allCellsX &&
           cellY >= 0 && cellY < allCellsY &&
           cellZ >= 0 && cellZ < allCellsZ;

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

    void CellIntersections (LinkedList<Intersection> result, Vector3d gridPos, Vector3d gridDir, Vector3d p0, Vector3d p1, CellRayIntersection cell)
    {
      void TriangleIntersection (Vector3d c1, Vector3d c2, Vector3d c3, Vector3d normal)
      {
        double t0 = Geometry.RayTriangleIntersection(
                        ref gridPos,
                        ref gridDir,
                        ref c1,
                        ref c2,
                        ref c3,
                        out var uv);

        if (t0 == double.NegativeInfinity)
          return;

        var intersection = new Intersection(this)
        {
          T = t0,
          Enter = true,
          Front = true,
          NormalLocal = normal,
          CoordLocal = p0 + t0 * p1,
        };

        // If the ray is not perpendicular on the face, then we are leaving
        var dot = Vector3d.Dot(intersection.NormalLocal, p1);
        if (dot > 0)
        {
          intersection.Enter = false;
          intersection.Front = false;
        };

        // TODO: Replace with ordering add?
        result.AddLast(intersection);
      }

      int cx = cell.Cx;
      int cy = cell.Cy;
      int cz = cell.Cz;
      double minT = cell.MinT;
      double maxT = cell.MaxT;
      gridPos -= new Vector3d(cx, cy, cz); // now gridPos is relative to cube of size 1; from (0,0,0) to (1,1,1)

      int corners = 0;
      double[,,] cube = new double[2,2,2];
      for (int i = 0; i <= 1; i++)
        for (int j = 0; j <= 1; j++)
          for (int k = 0; k <= 1; k++)
          {
            cube[i, j, k] = grid[cx + i, cy + j, cz + k] - 0.5;

            if (cube[i, j, k] > 0)
            {
              cube[i, j, k] = Math.Max(cube[i, j, k], 0.25);
              corners++;
            }
            else
              cube[i, j, k] = 0;
          }

      if (corners >= 8) // 8 corners -> do nothing
        return;
      if (corners == 5)
        InverseTriangle();
      if (corners == 6) // Only one line allowed
        InverseLines();
      if (corners > 4)
        InverseCorners();
      if (corners == 4)
        FourCornerMagnet();
      if (corners == 4)
        Four2x2Plane();
      if (corners == 4)
        FourSnake();
      if (corners >= 3)
        Triangles();
      if (corners >= 2)
        Lines();
      if (corners >= 1)
        Corners();

      // 4 points in a line
      void FourSnake ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
          {
            int k = 0;
            // "Snake" is always going to be as an order of 4 directions. We just need to find those directions.

            double r1 = cube[i,j,k];
            if (r1 == 0) // Check first corner
              continue;

            // Current corner
            int ci = i;
            int cj = j;
            int ck = k;

            // Opposite corners in each direction
            int i2 = (1-i);
            int j2 = (1-j);
            int k2 = (1-k);

            List<int> dirs = new List<int>();
            List<double> rs = new List<double>(){ r1 };
            // 1 = change in i
            // 2 = change in j
            // 3 = change in k

            for (int n = 0; n < 3; n++)
            {
              // We have to change i, j, or k
              if (i2 != -1 && cube[i2, cj, ck] > 0)
              {
                ci = i2;
                i2 = -1; // disable i2
                dirs.Add(1);
                rs.Add(cube[ci, cj, ck]);
                continue;
              }

              if (j2 != -1 && cube[ci, j2, ck] > 0)
              {
                cj = j2;
                j2 = -1; // disable i2
                dirs.Add(2);
                rs.Add(cube[ci, cj, ck]);
                continue;
              }

              if (k2 != -1 && cube[ci, cj, k2] > 0)
              {
                ck = k2;
                k2 = -1; // disable i2
                dirs.Add(3);
                rs.Add(cube[ci, cj, ck]);
                continue;
              }

              break;
            }

            // Check if we found the snake
            if (dirs.Count < 3)
              continue;

            Vector3d dirToVector (int direction)
            {
              switch (direction)
              {
                case 1:
                  return new Vector3d(1 - 2 * i, 0, 0);
                case 2:
                  return new Vector3d(0, 1 - 2 * j, 0);
                case 3:
                default:
                  return new Vector3d(0, 0, 1 - 2 * k);
              }
            }

            List<Vector3d> vertices = new List<Vector3d>();
            Vector3d start = new Vector3d(i,j,k);
            Vector3d dir1 = dirToVector(dirs[0]);
            Vector3d dir2 = dirToVector(dirs[1]);
            Vector3d dir3 = dirToVector(dirs[2]);
            Vector3d v0 = (start + dir3 * rs[0]);
            Vector3d v1 = (start + dir2 * rs[0]);
            Vector3d v2 = (start + dir1 + dir2 - dir1 * rs[2]);
            Vector3d v3 = (start + dir1 + dir2 + dir3 - dir1 * rs[3]);
            Vector3d v4 = (start + dir1 + dir2 + dir3 - dir2 * rs[3]);
            Vector3d v5 = (start + dir1 + dir2 * rs[1]);

            // Approximations
            Vector3d normal1 = dir2 + dir3;
            Vector3d normal2 = - dir1 - dir2;

            TriangleIntersection(v0, v1, v2, normal1); // It
            TriangleIntersection(v0, v2, v5, normal1); // It
            TriangleIntersection(v2, v5, v3, normal2); // It
            TriangleIntersection(v5, v3, v4, normal2); // It

            corners = 0;
          }
      }

      // Assumes there are no inverse lines or inverse triangles
      void InverseCorners ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
            for (int k = 0; k <= 1; k++)
            {
              double r = cube[i,j,k];
              double dirX = -2*i + 1; // Direction inside the cube
              double dirY = -2*j + 1;
              double dirZ = -2*k + 1;

              if (r <= 0)
              {
                Vector3d c1 = new Vector3d((1-i) - dirX * cube[(1-i),j,k], j                                     , k                                     );
                Vector3d c2 = new Vector3d(i                             , (1-j) + dirY * cube[i,(1-j),k]        , k                                     );
                Vector3d c3 = new Vector3d(i                             , j                                     , (1-k) + dirZ * cube[i,j,(1-k)]        );

                TriangleIntersection(c1, c2, c3, new Vector3d(-dirX, -dirY, -dirZ)); // It's just an approximation
              }
            }
        corners = 0;
        return;
      }

      void InverseLines ()
      {
        for (int dir = 0; dir < 3; dir++)
        {
          for (int i = 0; i <= 1; i++)
            for (int j = 0; j <= 1; j++)
            {
              // First corner coords
              int x1 = 0;
              int y1 = 0;
              int z1 = 0;

              // Second corner coords
              int x2 = 0;
              int y2 = 0;
              int z2 = 0;

              // Directions for both corners (for each direction)
              int dx1 = 0;
              int dy1 = 0;
              int dz1 = 0;
              int dx2 = 0;
              int dy2 = 0;
              int dz2 = 0;

              double dirI = (-2 * i + 1);
              double dirJ = (-2 * j + 1);

              Vector3d dir1 = Vector3d.Zero;
              Vector3d dir2 = Vector3d.Zero;

              switch (dir)
              {
                case 0: // Line along Z-axis
                {
                  x1 = x2 = i;
                  y1 = y2 = j;
                  z1 = 0;
                  z2 = 1;

                  dx1 = (int)dirI;
                  dy2 = (int)dirJ;

                  dir1 = Vector3d.UnitX * dirI;
                  dir2 = Vector3d.UnitY * dirJ;
                }
                break;
                case 1: // Line along Y-axis
                {
                  x1 = x2 = i;
                  y1 = 0;
                  y2 = 1;
                  z1 = z2 = j;

                  dx1 = (int)dirI;
                  dz2 = (int)dirJ;

                  dir1 = Vector3d.UnitX * dirI;
                  dir2 = Vector3d.UnitZ * dirJ;
                }
                break;
                case 2: // Line along X-Axis
                {
                  x1 = 0;
                  x2 = 1;
                  y1 = y2 = i;
                  z1 = z2 = j;

                  dy1 = (int)dirI;
                  dz2 = (int)dirJ;

                  dir1 = Vector3d.UnitY * dirI;
                  dir2 = Vector3d.UnitZ * dirJ;
                }
                break;
              }

              double r1 = cube[x1,y1,z1];
              double r2 = cube[x2,y2,z2];

              if (r1 <= 0 && r2 <= 0)
              {
                Vector3d c1 = new Vector3d(x1,y1,z1) + dir1 * (1 - cube[x1 + dx1, y1 + dy1, z1 + dz1]);
                Vector3d c2 = new Vector3d(x1,y1,z1) + dir2 * (1 - cube[x1 + dx2, y1 + dy2, z1 + dz2]);
                Vector3d c4 = new Vector3d(x2,y2,z2) + dir1 * (1 - cube[x2 + dx1, y2 + dy1, z2 + dz1]);
                Vector3d c3 = new Vector3d(x2,y2,z2) + dir2 * (1 - cube[x2 + dx2, y2 + dy2, z2 + dz2]);

                TriangleIntersection(c1, c2, c3, -dir1 - dir2);
                TriangleIntersection(c1, c3, c4, -dir1 - dir2);

                corners = 0;
                return;
                //// Handle the points
                //cube[x1, y1, z1] = 0;
                //cube[x2, y2, z2] = 0;
                //corners -= 2;
              }
            }
        }
      }

      // Inverse triangle (5 points)
      void InverseTriangle ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
            for (int k = 0; k <= 1; k++)
            {
              double r = cube[i,j,k];

              // Main corner is big
              if (r > 0)
                continue;

              // Opposite corners in each direction
              int i2 = (1-i);
              int j2 = (1-j);
              int k2 = (1-k);

              double ri = cube[i2,j,k];
              double rj = cube[i,j2,k];
              double rk = cube[i,j,k2];

              // Count if we have 2 neighbourhood points
              int count = (ri <= 0 ? 1 : 0) + (rj <= 0 ? 1 : 0) + (rk <= 0 ? 1 : 0);
              if (count < 2)
                continue;

              double dirI = -2*i + 1; // Direction from i2 to other side of the cube
              double dirJ = -2*j + 1;
              double dirK = -2*k + 1;

              Vector3d centerPoint;
              if (ri > 0)
                centerPoint = new Vector3d(i2 - dirI * ri, j, k);
              else if (rj > 0)
                centerPoint = new Vector3d(i, j2 - dirJ * rj, k);
              else
                centerPoint = new Vector3d(i, j, k2 - dirK * rk);

              List<Vector3d> vertices = new List<Vector3d>();

              if (ri <= 0)
              {
                if (rj > 0)
                {
                  vertices.Add(new Vector3d(i2, j2 - dirJ * cube[i2, j2, k], k));
                  vertices.Add(new Vector3d(i2, j, k2 - dirK * cube[i2, j, k2]));
                }
                else
                {
                  vertices.Add(new Vector3d(i2, j, k2 - dirK * cube[i2, j, k2]));
                  vertices.Add(new Vector3d(i2, j2 - dirJ * cube[i2, j2, k], k));
                }
              }
              if (rj <= 0)
              {
                if (rk > 0)
                {
                  vertices.Add(new Vector3d(i, j2, k2 - dirK * cube[i, j2, k2]));
                  vertices.Add(new Vector3d(i2 - dirI * cube[i2, j2, k], j2, k));
                }
                else
                {
                  vertices.Add(new Vector3d(i2 - dirI * cube[i2, j2, k], j2, k));
                  vertices.Add(new Vector3d(i, j2, k2 - dirK * cube[i, j2, k2]));
                }
              }
              if (rk <= 0)
              {
                if (rj > 0)
                {
                  vertices.Add(new Vector3d(i, j2 - dirJ * cube[i, j2, k2], k2));
                  vertices.Add(new Vector3d(i2 - dirI * cube[i2, j, k2], j, k2));
                }
                else
                {
                  vertices.Add(new Vector3d(i2 - dirI * cube[i2, j, k2], j, k2));
                  vertices.Add(new Vector3d(i, j2 - dirJ * cube[i, j2, k2], k2));
                }
              }

              Vector3d normal1;
              Vector3d normal2;
              normal2 = new Vector3d(-dirI, -dirJ, -dirK);
              if (ri > 0)
                normal1 = new Vector3d(-dirI, 0, 0);
              else if (rj > 0)
                normal1 = new Vector3d(0, -dirJ, 0);
              else
                normal1 = new Vector3d(0, 0, -dirK);

              TriangleIntersection(centerPoint, vertices[0], vertices[2], normal1);
              TriangleIntersection(vertices[0], vertices[1], vertices[2], normal2);
              TriangleIntersection(vertices[1], vertices[2], vertices[3], normal2);

              // handle
              corners = 0;
              return;
            }
      }

      // 4 points where 1 point is neighbour with all the others
      void FourCornerMagnet ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
            for (int k = 0; k <= 1; k++)
            {
              double r = cube[i,j,k];

              // Main corner is big
              if (r == 0)
                continue;

              // Opposite corners in each direction
              int i2 = (1-i);
              int j2 = (1-j);
              int k2 = (1-k);

              double ri = cube[i2,j,k];
              double rj = cube[i,j2,k];
              double rk = cube[i,j,k2];

              // Next corners are big
              if (ri == 0 || rj == 0 || rk == 0)
                continue;

              double dirI = -2*i + 1; // Direction from i2 to other side of the cube
              double dirJ = -2*j + 1;
              double dirK = -2*k + 1;

              // Names: c_<first direction fully>_<second direction by corner radius>
              Vector3d cij = new Vector3d(i2           , j + dirJ * ri, k            );
              Vector3d cik = new Vector3d(i2           , j            , k + dirK * ri);
              Vector3d cji = new Vector3d(i + dirI * rj, j2           , k            );
              Vector3d cjk = new Vector3d(i            , j2           , k + dirK * rj);
              Vector3d ckj = new Vector3d(i            , j + dirJ * rk, k2           );
              Vector3d cki = new Vector3d(i + dirI * rk, j            , k2           );

              Vector3d normal = new Vector3d(dirI, dirJ, dirK);

              TriangleIntersection(cij, cik, cji, normal);
              TriangleIntersection(cik, cji, cjk, normal);
              TriangleIntersection(cik, cjk, cki, normal);
              TriangleIntersection(cjk, cki, ckj, normal);

              // handle
              corners = 0;
              return;
            }
      }

      // quad from 4 points
      void Four2x2Plane ()
      {
        for (int face = 0; face < 6; face++)
        {
          // Variables made so that iterating i in 0..1 and j in 0..1 gives all corners of the face
          // Just set properly
          int bx = 0;
          int by = 0;
          int bz = 0;
          int xi = 0;
          int yi = 0;
          int yj = 0;
          int zj = 0;

          Vector3d normal = Vector3d.Zero;

          switch (face)
          {
            case 0: // z == 0
            {
              xi = 1;
              yj = 1;
              normal = new Vector3d(0, 0, 1);
            }
            break;
            case 1: // z == 1
            {
              bz = 1;
              xi = 1;
              yj = 1;
              normal = new Vector3d(0, 0, -1);
            }
            break;
            case 2: // x == 0
            {
              yi = 1;
              zj = 1;
              normal = new Vector3d(1, 0, 0);
            }
            break;
            case 3: // x == 1
            {
              bx = 1;
              yi = 1;
              zj = 1;
              normal = new Vector3d(-1, 0, 0);
            }
            break;
            case 4: // y == 0
            {
              xi = 1;
              zj = 1;
              normal = new Vector3d(0, 1, 0);
            }
            break;
            case 5: // y == 1
            {
              by = 1;
              xi = 1;
              zj = 1;
              normal = new Vector3d(0, -1, 0);
            }
            break;
          }

          // Count corners
          int c = 0;
          int ridx = 0;
          double[] rs = new double[4];

          for (int i = 0; i <= 1; i++)
            for (int j = 0; j <= 1; j++)
            {
              double r = cube[bx + xi * i, by + yi * i + yj * j, bz + zj * j];
              rs[ridx++] = r;
              if (r > 0)
                c++;
            }

          // If all corners are correct, make plane
          if (c == 4)
          {
            // ...(i,j)
            Vector3d c0 = new Vector3d(bx, by, bz) + normal * rs[0]; // (0,0)
            Vector3d c1 = new Vector3d(bx, by + yj, bz + zj) + normal * rs[1]; // (0,1)
            Vector3d c2 = new Vector3d(bx + xi, by + yi, bz) + normal * rs[2]; // (1,0)
            Vector3d c3 = new Vector3d(bx + xi, by + yi + yj, bz + zj) + normal * rs[3]; // (1,1)
            Vector3d cMid = new Vector3d(bx + 0.5 * xi, by + 0.5 * yi + 0.5 * yj, bz + 0.5 * zj) + normal * (rs[0] + rs[1] + rs[2] + rs[3]) * 0.25; // (1,1)

            TriangleIntersection(c0, c1, cMid, normal);
            TriangleIntersection(c1, c3, cMid, normal);
            TriangleIntersection(c3, c2, cMid, normal);
            TriangleIntersection(c2, c0, cMid, normal);

            // Handle the points
            //cube[bx, by, bz] = 0;
            //cube[bx, by + yj, bz + zj] = 0;
            //cube[bx + xi, by + yi, bz] = 0;
            //cube[bx + xi, by + yi + yj, bz + zj] = 0;
            //corners -= 4;
            corners = 0;
            return;
          }
        }
      }

      void Triangles ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
            for (int k = 0; k <= 1; k++)
            {
              double r = cube[i,j,k];

              // Main corner is big
              if (r == 0)
                continue;

              // Opposite corners in each direction
              int i2 = (1-i);
              int j2 = (1-j);
              int k2 = (1-k);

              double ri = cube[i2,j,k];
              double rj = cube[i,j2,k];
              double rk = cube[i,j,k2];

              // Count if we have 2 neighbourhood points
              int count = (ri > 0 ? 1 : 0) + (rj > 0 ? 1 : 0) + (rk > 0 ? 1 : 0);
              if (count < 2)
                continue;

              double dirI = -2*i + 1; // Direction from i2 to other side of the cube
              double dirJ = -2*j + 1;
              double dirK = -2*k + 1;

              Vector3d centerPoint;
              if (ri <= 0)
                centerPoint = new Vector3d(i + dirI * r, j, k);
              else if (rj <= 0)
                centerPoint = new Vector3d(i, j + dirJ * r, k);
              else
                centerPoint = new Vector3d(i, j, k + dirK * r);

              List<Vector3d> vertices = new List<Vector3d>();

              if (ri > 0)
              {
                if (rj <= 0)
                {
                  vertices.Add(new Vector3d(i2, j + dirJ * ri, k));
                  vertices.Add(new Vector3d(i2, j, k + dirK * ri));
                }
                else
                {
                  vertices.Add(new Vector3d(i2, j, k + dirK * ri));
                  vertices.Add(new Vector3d(i2, j + dirJ * ri, k));
                }
              }
              if (rj > 0)
              {
                if (rk <= 0)
                {
                  vertices.Add(new Vector3d(i, j2, k + dirK * rj));
                  vertices.Add(new Vector3d(i + dirI * rj, j2, k));
                }
                else
                {
                  vertices.Add(new Vector3d(i + dirI * rj, j2, k));
                  vertices.Add(new Vector3d(i, j2, k + dirK * rj));
                }
              }
              if (rk > 0)
              {
                if (rj <= 0)
                {
                  vertices.Add(new Vector3d(i, j + dirJ * rk, k2));
                  vertices.Add(new Vector3d(i + dirI * rk, j, k2));
                }
                else
                {
                  vertices.Add(new Vector3d(i + dirI * rk, j, k2));
                  vertices.Add(new Vector3d(i, j + dirJ * rk, k2));
                }
              }

              Vector3d normal1;
              Vector3d normal2;
              if (ri <= 0)
              {
                normal1 = new Vector3d(dirI, 0, 0);
                normal2 = new Vector3d(dirI, dirJ, dirK);
              }
              else if (rj <= 0)
              {
                normal1 = new Vector3d(0, dirJ, 0);
                normal2 = new Vector3d(dirI, dirJ, dirK);
              }
              else
              {
                normal1 = new Vector3d(0, 0, dirK);
                normal2 = new Vector3d(dirI, dirJ, dirK);
              }

              TriangleIntersection(centerPoint, vertices[0], vertices[2], normal1);
              TriangleIntersection(vertices[0], vertices[1], vertices[2], normal2);
              TriangleIntersection(vertices[1], vertices[2], vertices[3], normal2);

              // handle
              corners -= 3;
              if (corners < 3)
                return;
            }
      }

      void Lines ()
      {
        for (int dir = 0; dir < 3; dir++)
        {
          for (int i = 0; i <= 1; i++)
            for (int j = 0; j <= 1; j++)
            {
              // First corner coords
              int x1 = 0;
              int y1 = 0;
              int z1 = 0;

              // Second corner coords
              int x2 = 0;
              int y2 = 0;
              int z2 = 0;

              double dirI = (-2 * i + 1);
              double dirJ = (-2 * j + 1);

              Vector3d dir1 = Vector3d.Zero;
              Vector3d dir2 = Vector3d.Zero;

              switch (dir)
              {
                case 0:
                {
                  x1 = x2 = i;
                  y1 = y2 = j;
                  z1 = 0;
                  z2 = 1;

                  dir1 = Vector3d.UnitX * dirI;
                  dir2 = Vector3d.UnitY * dirJ;
                }
                break;
                case 1:
                {
                  x1 = x2 = i;
                  y1 = 0;
                  y2 = 1;
                  z1 = z2 = j;

                  dir1 = Vector3d.UnitX * dirI;
                  dir2 = Vector3d.UnitZ * dirJ;
                }
                break;
                case 2:
                {
                  x1 = 0;
                  x2 = 1;
                  y1 = y2 = i;
                  z1 = z2 = j;

                  dir1 = Vector3d.UnitY * dirI;
                  dir2 = Vector3d.UnitZ * dirJ;
                }
                break;
              }

              double r1 = cube[x1,y1,z1];
              double r2 = cube[x2,y2,z2];

              if (r1 > 0 && r2 > 0)
              {
                Vector3d c1 = new Vector3d(x1,y1,z1) + dir1 * r1;
                Vector3d c2 = new Vector3d(x1,y1,z1) + dir2 * r1;
                Vector3d c4 = new Vector3d(x2,y2,z2) + dir1 * r2;
                Vector3d c3 = new Vector3d(x2,y2,z2) + dir2 * r2;

                TriangleIntersection(c1, c2, c3, dir1 + dir2);
                TriangleIntersection(c1, c3, c4, dir1 + dir2);

                // Handle the points
                cube[x1, y1, z1] = 0;
                cube[x2, y2, z2] = 0;
                corners -= 2;
                if (corners < 2)
                  return;
              }
            }
        }
      }

      void Corners ()
      {
        for (int i = 0; i <= 1; i++)
          for (int j = 0; j <= 1; j++)
            for (int k = 0; k <= 1; k++)
            {
              double r = cube[i,j,k];
              double dirX = -2*i + 1; // Direction inside the cube
              double dirY = -2*j + 1;
              double dirZ = -2*k + 1;

              if (r > 0)
              {
                Vector3d c1 = new Vector3d(i + dirX * r, j           , k);
                Vector3d c2 = new Vector3d(i           , j + dirY * r, k);
                Vector3d c3 = new Vector3d(i           , j           , k + dirZ * r);

                TriangleIntersection(c1, c2, c3, new Vector3d(dirX, dirY, dirZ));
              }
            }
      }
    }

    double t = -1;
    public override LinkedList<Intersection> Intersect (Vector3d p0, Vector3d p1)
    {
      if (Time > t)
      {
        Debug.WriteLine(Time);
        t = Time;
      }
      LinkedList<Intersection> result = new LinkedList<Intersection>();

      Vector3d gridPos = new Vector3d((p0.X + 0.5) * allCellsX, (p0.Y + 0.5) * allCellsY, (p0.Z + 0.5) * allCellsZ);
      Vector3d gridDir = new Vector3d(p1.X * allCellsX, p1.Y * allCellsY, p1.Z * allCellsZ);

      foreach (var cri in GridRayIntersectingCells(gridPos, gridDir))
        CellIntersections(result, gridPos, gridDir, p0, p1, cri);

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
