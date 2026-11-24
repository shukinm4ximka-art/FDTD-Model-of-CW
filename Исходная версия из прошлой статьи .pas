program FDTD3D_Animated_PEC_Energy;

uses
  GraphABC, System;

type
  Real = double;

const
  Nx = 100;
  Ny = 100;
  Nz = 100;

  dx = 1.0e-3;
  dy = 1.0e-3;
  dz = 1.0e-3;
  c  = 2.998e8;
  eps0 = 8.854e-12;
  mu0  = 4.0 * Pi * 1.0e-7;

  steps = 300;
  snap_every = 5;

  t0 = 30.0;
  sigma_t = 10.0;

  J0 = 1.0e6;

  // ---------- Раскладка окна ----------
  WinWidth  = 1240;
  WinHeight = 660;
  AnimX     = 10;
  AnimY     = 40;
  AnimSize  = 600;
  PlotX     = 630;
  PlotY     = 40;
  PlotW     = 590;
  PlotH     = 600;

var
  dt: Real;
  Ex, Ey, Ez, Hx, Hy, Hz: array of Real;
  eps_r, mu_r: array of Real;
  src_i, src_j, src_k, slice_z: integer;
  W, H: integer;
  E_hist: array of Real;   // история полной энергии

// ---------- Индексаторы Yee-сетки ----------
function idxEx(i, j, k: integer): integer;
begin Result := (i * (Ny - 1) + j) * (Nz - 1) + k; end;
function idxEy(i, j, k: integer): integer;
begin Result := (i * Ny + j) * (Nz - 1) + k; end;
function idxEz(i, j, k: integer): integer;
begin Result := (i * (Ny - 1) + j) * Nz + k; end;
function idxHx(i, j, k: integer): integer;
begin Result := (i * Ny + j) * Nz + k; end;
function idxHy(i, j, k: integer): integer;
begin Result := (i * (Ny - 1) + j) * Nz + k; end;
function idxHz(i, j, k: integer): integer;
begin Result := (i * Ny + j) * (Nz - 1) + k; end;
function idxEps(i, j, k: integer): integer;
begin Result := (i * Ny + j) * Nz + k; end;

// ---------- Источник ----------
function source(n: integer): Real;
begin
  Result := Exp(-Sqr((n - t0) / sigma_t));
end;

// ---------- PEC: E_tang = 0 на нижних гранях ----------
procedure ApplyPEC;
var i, j, k: integer;
begin
  for j := 0 to Ny - 1 do
    for k := 0 to Nz - 2 do
      Ey[idxEy(0, j, k)] := 0;
  for j := 0 to Ny - 2 do
    for k := 0 to Nz - 1 do
      Ez[idxEz(0, j, k)] := 0;

  for i := 0 to Nx - 1 do
    for k := 0 to Nz - 2 do
      Ex[idxEx(i, 0, k)] := 0;
  for i := 0 to Nx - 2 do
    for k := 0 to Nz - 1 do
      Ez[idxEz(i, 0, k)] := 0;

  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 2 do
      Ex[idxEx(i, j, 0)] := 0;
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 1 do
      Ey[idxEy(i, j, 0)] := 0;
end;

// ---------- Полная энергия ----------
function TotalEnergy: Real;
var i, j, k: integer; sum: Real;
begin
  sum := 0;
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 2 do
        sum += eps0 * Sqr(Ex[idxEx(i,j,k)]);
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 2 do
        sum += eps0 * Sqr(Ey[idxEy(i,j,k)]);
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 1 do
        sum += eps0 * Sqr(Ez[idxEz(i,j,k)]);
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 1 do
        sum += mu0 * Sqr(Hx[idxHx(i,j,k)]);
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 1 do
        sum += mu0 * Sqr(Hy[idxHy(i,j,k)]);
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 2 do
        sum += mu0 * Sqr(Hz[idxHz(i,j,k)]);
  Result := 0.5 * sum * dx * dy * dz;
end;

// ---------- max |div E| ----------
function MaxDivergenceE: Real;
var i, j, k: integer; d, maxd: Real;
begin
  maxd := 0;
  for i := 1 to Nx - 2 do
    for j := 1 to Ny - 2 do
      for k := 1 to Nz - 2 do
      begin
        d := (Ex[idxEx(i+1,j,k)] - Ex[idxEx(i,j,k)]) / dx
           + (Ey[idxEy(i,j+1,k)] - Ey[idxEy(i,j,k)]) / dy
           + (Ez[idxEz(i,j,k+1)] - Ez[idxEz(i,j,k)]) / dz;
        if Abs(d) > maxd then maxd := Abs(d);
      end;
  Result := maxd;
end;

// ---------- Срез Ez (рисует только свою панель) ----------
procedure DrawSlice(frameArr: array of Real; x0, y0, sz: integer);
var
  i, j, cell: integer;
  val, norm: Real;
  v: byte;
  col: Color;
begin
  norm := 0.0;
  for i := 0 to W * H - 1 do
    if Abs(frameArr[i]) > norm then norm := Abs(frameArr[i]);
  if norm = 0 then norm := 1.0;

  cell := sz div W;
  if cell < 1 then cell := 1;

  for i := 0 to W - 1 do
    for j := 0 to H - 1 do
    begin
       val := frameArr[i * H + j]/norm;
       if val < 0
           then col := RGB(255, 0, 255 + Round(255 * val))
       else
           col := RGB(255, 255 - Round(255 * val), 0);
      Brush.Color := col;
      FillRect(x0 + i * cell, y0 + j * cell,
               x0 + (i + 1) * cell, y0 + (j + 1) * cell);
    end;

  Font.Color := clBlack;
  Font.Size := 11;
  TextOut(x0 + 5, y0 - 22, 'Срез Ez, z = ' + IntToStr(slice_z));
end;

// ---------- График энергии W(t) ----------
procedure DrawEnergyPlot(E_hist: array of Real; nCur, x0, y0, w, h: integer);
var
  i, px, py: integer;
  Emax, frac: Real;
  marginL, marginB, marginT, marginR, plotW, plotH: integer;
begin
  marginL := 60;
  marginB := 40;
  marginT := 25;
  marginR := 20;
  plotW := w - marginL - marginR;
  plotH := h - marginT - marginB;

  Emax := 1e-30;
  for i := 0 to nCur do
    if E_hist[i] > Emax then Emax := E_hist[i];
  if Emax <= 0 then Emax := 1;

  // Оси
  Pen.Width := 1;
  Pen.Color := clBlack;
  Line(x0 + marginL, y0 + marginT, x0 + marginL, y0 + marginT + plotH);
  Line(x0 + marginL, y0 + marginT + plotH,
       x0 + marginL + plotW, y0 + marginT + plotH);

  Font.Color := clBlack;
  Font.Size := 9;

  // Y-метки: 0%, 25%, 50%, 75%, 100%
  for i := 0 to 4 do
  begin
    frac := i / 4;
    py := y0 + marginT + plotH - Round(plotH * frac);
    Line(x0 + marginL - 5, py, x0 + marginL, py);
    TextOut(x0 + 10, py - 7, IntToStr(Round(frac * 100)) + '%');
  end;

  // X-метки
  for i := 0 to 4 do
  begin
    px := x0 + marginL + Round(plotW * i / 4);
    Line(px, y0 + marginT + plotH, px, y0 + marginT + plotH + 5);
    TextOut(px - 12, y0 + marginT + plotH + 8,
            IntToStr(Round((steps - 1) * i / 4)));
  end;

  // Кривая
  Pen.Width := 2;
  Pen.Color := clRed;
  for i := 0 to nCur do
  begin
    px := x0 + marginL + Round(plotW * i / (steps - 1));
    py := y0 + marginT + plotH - Round(plotH * E_hist[i] / Emax);
    if i = 0 then MoveTo(px, py) else LineTo(px, py);
  end;

  // «Бегунок» текущего шага
  if nCur > 0 then
  begin
    Pen.Width := 1;
    Pen.Color := clBlue;
    px := x0 + marginL + Round(plotW * nCur / (steps - 1));
    Line(px, y0 + marginT, px, y0 + marginT + plotH);
  end;

  // Заголовок и абсолютный максимум
  Font.Color := clBlack;
  Font.Size := 11;
  TextOut(x0 + marginL, y0 + 2, 'Полная энергия W(t)');
  Font.Size := 9;
  TextOut(x0 + marginL + 5, y0 + marginT + 5, Format('W_max = {0:E2}', Emax));
  TextOut(x0 + marginL + plotW div 2 - 25, y0 + h - 15, 'номер шага');
end;

// ======================== ОСНОВНАЯ ПРОГРАММА ========================
begin
  SetWindowSize(WinWidth, WinHeight);
  Window.Title := 'FDTD 3D PEC: срез Ez и энергия';

  dt := 0.95 * dx / (Sqrt(3.0) * c);
  WriteLn('CFL: dt = ', dt, '  c*dt/dx = ', c*dt/dx,
          '  (предел 1/sqrt(3) = 0.5774)');

  src_i := Nx div 2;
  src_j := Ny div 2;
  src_k := Nz div 2;
  slice_z := Nz div 2;

  SetLength(Ex, Nx * (Ny - 1) * (Nz - 1));
  SetLength(Ey, (Nx - 1) * Ny * (Nz - 1));
  SetLength(Ez, (Nx - 1) * (Ny - 1) * Nz);
  SetLength(Hx, (Nx - 1) * Ny * Nz);
  SetLength(Hy, Nx * (Ny - 1) * Nz);
  SetLength(Hz, Nx * Ny * (Nz - 1));
  SetLength(eps_r, Nx * Ny * Nz);
  SetLength(mu_r,  Nx * Ny * Nz);
  SetLength(E_hist, steps);

  for var i := 0 to High(eps_r) do
  begin
    eps_r[i] := 1.0;
    mu_r[i]  := 1.0;
  end;

  W := Nx - 1;
  H := Ny - 1;
  var frameArr := new Real[W * H];

  // ================= Временной цикл =================
  for var n := 0 to steps - 1 do
  begin
            
        // ---------- Обновление H (ИСПРАВЛЕНО: запись в собственный слот) ----------
    for var i := 0 to Nx - 2 do
      for var j := 1 to Ny - 2 do
        for var k := 1 to Nz - 2 do
        begin
          var dEz_dy := (Ez[idxEz(i, j, k)] - Ez[idxEz(i, j - 1, k)]) / dy;
          var dEy_dz := (Ey[idxEy(i, j, k)] - Ey[idxEy(i, j, k - 1)]) / dz;
          var mu := mu0 * mu_r[idxEps(i, j, k)];
          Hx[idxHx(i, j, k)] += -(dEz_dy - dEy_dz) / mu * dt;
        end;

    for var i := 1 to Nx - 2 do
      for var j := 0 to Ny - 2 do
        for var k := 1 to Nz - 2 do
        begin
          var dEx_dz := (Ex[idxEx(i, j, k)] - Ex[idxEx(i, j, k - 1)]) / dz;
          var dEz_dx := (Ez[idxEz(i, j, k)] - Ez[idxEz(i - 1, j, k)]) / dx;
          var mu := mu0 * mu_r[idxEps(i, j, k)];
          Hy[idxHy(i, j, k)] += -(dEx_dz - dEz_dx) / mu * dt;
        end;

    for var i := 1 to Nx - 2 do
     for var j := 1 to Ny - 2 do
        for var k := 0 to Nz - 2 do
        begin
          var dEy_dx := (Ey[idxEy(i, j, k)] - Ey[idxEy(i - 1, j, k)]) / dx;
          var dEx_dy := (Ex[idxEx(i, j, k)] - Ex[idxEx(i, j - 1, k)]) / dy;
          var mu := mu0 * mu_r[idxEps(i, j, k)];
          Hz[idxHz(i, j, k)] += -(dEy_dx - dEx_dy) / mu * dt;
        end;

    // ---------- Обновление E ----------
    for var i := 0 to Nx - 1 do
      for var j := 0 to Ny - 2 do
        for var k := 0 to Nz - 2 do
        begin
          var dHz_dy := (Hz[idxHz(i, j + 1, k)] - Hz[idxHz(i, j, k)]) / dy;
          var dHy_dz := (Hy[idxHy(i, j, k + 1)] - Hy[idxHy(i, j, k)]) / dz;
          var eps := eps0 * eps_r[idxEps(i, j, k)];
          Ex[idxEx(i, j, k)] += (dHz_dy - dHy_dz) / eps * dt;
          if (i = src_i) and (j = src_j - 1) and (k = src_k - 1) then
            Ex[idxEx(i, j, k)] += dt * J0 * source(n) / eps;
        end;

    for var i := 0 to Nx - 2 do
      for var j := 0 to Ny - 1 do
        for var k := 0 to Nz - 2 do
        begin
          var dHx_dz := (Hx[idxHx(i, j, k + 1)] - Hx[idxHx(i, j, k)]) / dz;
          var dHz_dx := (Hz[idxHz(i + 1, j, k)] - Hz[idxHz(i, j, k)]) / dx;
          var eps := eps0 * eps_r[idxEps(i, j, k)];
          Ey[idxEy(i, j, k)] += (dHx_dz - dHz_dx) / eps * dt;
        end;

    for var i := 0 to Nx - 2 do
      for var j := 0 to Ny - 2 do
        for var k := 0 to Nz - 1 do
        begin
          var dHy_dx := (Hy[idxHy(i + 1, j, k)] - Hy[idxHy(i, j, k)]) / dx;
          var dHx_dy := (Hx[idxHx(i, j + 1, k)] - Hx[idxHx(i, j, k)]) / dy;
          var eps := eps0 * eps_r[idxEps(i, j, k)];
          Ez[idxEz(i, j, k)] += (dHy_dx - dHx_dy) / eps * dt;
        end;

    ApplyPEC;

    // ---------- Энергия на каждом шаге ----------
    E_hist[n] := TotalEnergy;

    // ---------- Диагностика и отрисовка ----------
    if (n mod snap_every = 0) or (n = steps - 1) then
    begin
      WriteLn('n=', n:4,
              '  W=', E_hist[n]:0:6,
              '  max|divE|=', MaxDivergenceE:0:6);

      for var i := 0 to W - 1 do
        for var j := 0 to H - 1 do
          frameArr[i * H + j] := Ez[idxEz(i, j, slice_z)];

      LockDrawing;
      ClearWindow(clWhite);
      DrawSlice(frameArr, AnimX, AnimY, AnimSize);
      DrawEnergyPlot(E_hist, n, PlotX, PlotY, PlotW, PlotH);
      Redraw;
      SaveWindow('frame_' + IntToStr(n div snap_every) + '.png');
      Sleep(10);
    end;
  end;

  WriteLn('Расчёт завершён.');
  WriteLn(Ez[idxEz(src_i - 10, src_j, slice_z)],
        Ez[idxEz(src_i + 10, src_j, slice_z)]);

end.