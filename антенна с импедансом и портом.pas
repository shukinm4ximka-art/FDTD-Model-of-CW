program FDTD3D_Dipole_CW;

uses
  GraphABC, System;

type
  Real = double;

const
  Nx = 30;
  Ny = 30;
  Nz = 30;

  dx = 1.0e-3;
  dy = 1.0e-3;
  dz = 1.0e-3;
  c  = 2.998e8;
  eps0 = 8.854e-12;
  mu0  = 4.0 * Pi * 1.0e-7;

  // ---------- ДОБАВЛЕНО: геометрия антенны ----------
  // Диполь вдоль x: плечи по ArmLen ячеек + зазор 1 ячейка в src_i.
  // Полная длина ~11 мм -> резонанс f0 = c/2L ~ 15.2 ГГц, lambda ~ 19.7 мм.
  // ДОБАВЛЕНО следствие: PeriodSteps = 36 -> ~20 ячеек на длину волны
  // (было 16 и 8.8 ячеек — претензия "грубая сетка" закрыта выбором частоты).
  ArmLen      = 5;      // плечо диполя, ячеек
  V0          = 10.0;   // амплитуда питающего напряжения, В

  // ---------- ИЗМЕНЕНО: параметры несущей под диполь ----------
  PeriodSteps = 36;     // шагов на период несущей (lambda ~ 20 ячеек)
  KeyRamp     = 200;    // плавный фронт ключа, шагов (~5.5 периода несущей)

  // ---------- морзянка ----------
  Msg = 'MAX2';         // сообщение любое: steps и snap_every пересчитаются сами
  unitSteps = 1000;     // длительность точки, шагов
  maxP = 3000;          // запас под расписание

  // ---------- длительность прогона (вычисляется из Msg) ----------
  TAIL         = 20000; // хвост: ждём честный спад энергии после конца передачи
  FRAME_TARGET = 60;    // желаемое число кадров анимации
  SNAP_ROUND   = 100;   // snap округляется вверх до этой сотни (ровные оси)

  // ---------- ИЗМЕНЕНО: губка сильнее — закрывает претензию о резонаторе ----------
  AbsDepth = 8;         // толщина слоя, ячеек (было 4)
  AbsMax   = 0.5;       // максимум потерь на внешней грани (было 0.2)

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

  // вычисляются из Msg, а не задаются константой
  steps, snap_every: integer;

  // ---------- CW-расписание: интервалы включения ключа ----------
  onStart, onEnd: array[1..maxP] of integer;
  nP: integer;

  dampX: array[0..Nx-1] of Real;
  dampY: array[0..Ny-1] of Real;
  dampZ: array[0..Nz-1] of Real;

  // ---------- ДОБАВЛЕНО: диагностика порта (входной импеданс) ----------
  pr: text;   // пробник
  gr: text;   // gap: V(t), I(t) в зазоре
  prevGap: Real;
  ZreV, ZimV, ZreI, ZimI: Real;   // накопители ДПФ V и I на несущей
  w0: Real;                       // omega несущей, рад/шаг

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

// ---------- Таблица морзянки ----------
// ЕДИНЫЙ ИСТОЧНИК ИСТИНЫ: и построитель расписания, и MsgUnits
// считают тайминги только через эту таблицу и unitSteps.
function MorseCode(ch: char): string;
begin
  case ch of
    'A': Result := '.-';   'B': Result := '-...'; 'C': Result := '-.-.';
    'D': Result := '-..';  'E': Result := '.';    'F': Result := '..-.';
    'G': Result := '--.';  'H': Result := '....'; 'I': Result := '..';
    'J': Result := '.---'; 'K': Result := '-.-';  'L': Result := '.-..';
    'M': Result := '--';   'N': Result := '-.';   'O': Result := '---';
    'P': Result := '.--.'; 'Q': Result := '--.-'; 'R': Result := '.-.';
    'S': Result := '...';  'T': Result := '-';    'U': Result := '..-';
    'V': Result := '...-'; 'W': Result := '.--';  'X': Result := '-..-';
    'Y': Result := '-.--'; 'Z': Result := '--..';
    '0': Result := '-----'; '1': Result := '.----';
    '2': Result := '..---'; '3': Result := '...--';
    '4': Result := '....-'; '5': Result := '.....';
    '6': Result := '-....'; '7': Result := '--...';
    '8': Result := '---..'; '9': Result := '----.';
    '.': Result := '.-.-.-';
    else Result := '';
  end;
end;

// ---------- Длительность сообщения в юнитах ----------
// Зеркало построителя расписания: элемент вносит (длина + 1 пауза),
// после знака +2 (итого 3 юнита между знаками), пробел = 5 юнитов.
function MsgUnits(const s: string): integer;
var
  i, k: integer;
  code: string;
begin
  Result := 0;
  for i := 1 to Length(s) do
  begin
    if s[i] = ' ' then
      Result += 5
    else
    begin
      code := MorseCode(s[i]);
      if code <> '' then
      begin
        for k := 1 to Length(code) do
          if code[k] = '-' then Result += 3 + 1 else Result += 1 + 1;
        Result += 2;   // пауза между знаками
      end;
    end;
  end;
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

// ---------- Инициализация профиля губки ----------
procedure InitSponge;
var i: integer; dd: Real;
begin
  for i := 0 to Nx - 1 do dampX[i] := 1.0;
  for i := 0 to Ny - 1 do dampY[i] := 1.0;
  for i := 0 to Nz - 1 do dampZ[i] := 1.0;
  for i := 0 to AbsDepth - 1 do
  begin
    dd := AbsMax * Sqr((AbsDepth - i) / AbsDepth);   // квадратичный профиль
    dampX[i] := 1 - dd;  dampX[Nx - 1 - i] := 1 - dd;
    dampY[i] := 1 - dd;  dampY[Ny - 1 - i] := 1 - dd;
    dampZ[i] := 1 - dd;  dampZ[Nz - 1 - i] := 1 - dd;
  end;
end;

// ---------- Применение губки (после PEC, до записи пробника) ----------
procedure ApplySponge;
var i, j, k: integer;
begin
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 2 do
        Ex[idxEx(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 2 do
        Ey[idxEy(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 1 do
        Ez[idxEz(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
  for i := 0 to Nx - 2 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 1 do
        Hx[idxHx(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 2 do
      for k := 0 to Nz - 1 do
        Hy[idxHy(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
  for i := 0 to Nx - 1 do
    for j := 0 to Ny - 1 do
      for k := 0 to Nz - 2 do
        Hz[idxHz(i, j, k)] *= dampX[i] * dampY[j] * dampZ[k];
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

// ---------- Срез Ez ----------
procedure DrawSlice(frameArr: array of Real; x0, y0, sz: integer);
var
  i, j, cell: integer;
  val, norm: Real;
  col: Color;
begin
  norm := 0.0;
  for i := 0 to W * H - 1 do
    if Abs(frameArr[i]) > norm then norm := Abs(frameArr[i]);
  if norm = 0 then norm := 1.0;
  norm := norm / 2;    // экспериментальное усиление;необходимо /20

  cell := sz div W;
  if cell < 1 then cell := 1;

  for i := 0 to W - 1 do
    for j := 0 to H - 1 do
    begin
      val := frameArr[i * H + j] / norm;
      if val < 0
          then col := RGB(255, 0, 255 + Round(255 * val))
          else col := RGB(255, 255 - Round(255 * val), 0);
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

  Pen.Width := 1;
  Pen.Color := clBlack;
  Line(x0 + marginL, y0 + marginT, x0 + marginL, y0 + marginT + plotH);
  Line(x0 + marginL, y0 + marginT + plotH,
       x0 + marginL + plotW, y0 + marginT + plotH);

  Font.Color := clBlack;
  Font.Size := 9;

  for i := 0 to 4 do
  begin
    frac := i / 4;
    py := y0 + marginT + plotH - Round(plotH * frac);
    Line(x0 + marginL - 5, py, x0 + marginL, py);
    TextOut(x0 + 10, py - 7, IntToStr(Round(frac * 100)) + '%');
  end;

  for i := 0 to 4 do
  begin
    px := x0 + marginL + Round(plotW * i / 4);
    Line(px, y0 + marginT + plotH, px, y0 + marginT + plotH + 5);
    TextOut(px - 12, y0 + marginT + plotH + 8,
            IntToStr(Round((steps - 1) * i / 4)));
  end;

  Pen.Width := 2;
  Pen.Color := clRed;
  for i := 0 to nCur do
  begin
    px := x0 + marginL + Round(plotW * i / (steps - 1));
    py := y0 + marginT + plotH - Round(plotH * E_hist[i] / Emax);
    if i = 0 then MoveTo(px, py) else LineTo(px, py);
  end;

  if nCur > 0 then
  begin
    Pen.Width := 1;
    Pen.Color := clBlue;
    px := x0 + marginL + Round(plotW * nCur / (steps - 1));
    Line(px, y0 + marginT, px, y0 + marginT + plotH);
  end;

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
  Window.Title := 'FDTD 3D: диполь-антенна, морзянка «' + Msg + '» (CW)';

  dt := 0.95 * dx / (Sqrt(3.0) * c);
  WriteLn('CFL: dt = ', dt, '  c*dt/dx = ', c*dt/dx,
          '  (предел 1/sqrt(3) = 0.5774)');
  WriteLn('Диполь: плечи ', ArmLen, ' ячеек, длина ~', 2*ArmLen + 1, ' мм;',
          ' f0 = ', (c / (2 * (2*ArmLen + 1) * dx)) / 1.0e9:0:2, ' ГГц,',
          ' ячеек на lambda = ', (PeriodSteps * 0.95 / Sqrt(3.0)):0:1);

  // ---------- steps и snap_every вычисляются из сообщения ----------
  steps := MsgUnits(Msg) * unitSteps + TAIL;
  snap_every := (steps + FRAME_TARGET - 1) div FRAME_TARGET;
  snap_every := ((snap_every + SNAP_ROUND - 1) div SNAP_ROUND) * SNAP_ROUND;
  if snap_every < SNAP_ROUND then snap_every := SNAP_ROUND;
  WriteLn('Сообщение «', Msg, '»: ', MsgUnits(Msg), ' юнитов -> steps = ', steps,
          ', snap_every = ', snap_every, ' (~', steps div snap_every, ' кадров)');

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

  // ---------- CW-расписание: интервалы включения ключа ----------
  nP := 0;
  var curPos := 0;
  for var k := 1 to Length(Msg) do
  begin
    if Msg[k] = ' ' then begin curPos := curPos + 5 * unitSteps; continue; end;
    var code := MorseCode(Msg[k]);
    for var s := 1 to Length(code) do
    begin
      var el := 1;
      if code[s] = '-' then el := 3;
      Inc(nP);
      if nP > maxP then
      begin
        WriteLn('ОШИБКА: расписание не влезло в onStart/onEnd, maxP = ', maxP);
        Halt;
      end;
      onStart[nP] := curPos;
      onEnd[nP]   := curPos + el * unitSteps - 1;
      curPos := curPos + (el + 1) * unitSteps;   // элемент + пауза внутри знака
    end;
    curPos := curPos + 2 * unitSteps;            // пауза между знаками (итого 3 юнита)
  end;
  WriteLn('CW-посылка: ', nP, ' интервалов, конец на шаге ', curPos,
          ', далее хвост ', TAIL, ' шагов');

  // ---------- Файлы: пробник и порт ----------
  Assign(pr, 'probe_cw_антенна.txt');
  Rewrite(pr);
  Assign(gr, 'gap.txt');
  Rewrite(gr);

  prevGap := 0.0;
  ZreV := 0; ZimV := 0; ZreI := 0; ZimI := 0;
  w0 := 2 * Pi / PeriodSteps;   // omega несущей, рад/шаг

  InitSponge;

  // ================= Временной цикл =================
  for var n := 0 to steps - 1 do
  begin

    // ---------- Обновление H ----------
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

          // ---------- ДОБАВЛЕНО: диполь — плечи PEC + зазор с источником ----------
          // Тонкий провод вдоль оси x: обнуляем продольную E на металле
          // (касательная E на поверхности проводника). Зазор — источник
          // напряжения: жёсткий при включённом ключе, свободный при выключенном
          // (иначе принудительный ноль в паузе "замкнул" бы зазор и убил звон).
          if (k = src_k - 1) and (j = src_j - 1) then
          begin
            for var m := 1 to ArmLen do
            begin
              if i = src_i - m then Ex[idxEx(i, j, k)] := 0;
              if i = src_i + m then Ex[idxEx(i, j, k)] := 0;
            end;
            if i = src_i then
            begin
              var sVal := 0.0;
              for var p := 1 to nP do
                if (n >= onStart[p]) and (n <= onEnd[p]) then
                begin
                  // огибающая ключа: линейный фронт и спад
                  var env := 1.0;
                  if n - onStart[p] < KeyRamp then
                    env := (n - onStart[p]) / KeyRamp;
                  if onEnd[p] - n < KeyRamp then
                  begin
                    var e2 := (onEnd[p] - n) / KeyRamp;
                    if e2 < env then env := e2;
                  end;
                  sVal += env * Sin(2 * Pi * (n - onStart[p]) / PeriodSteps);
                end;
              // жёсткий источник напряжения: E = V/dx (нулевой выходной импеданс).
              // Если sVal = 0 (пауза) — зазор "замкнут накоротко", это честная
              // модель замкнутого ключа: звон антенны продолжается, источник
              // просто шунтирует зазор.
              Ex[idxEx(i, j, k)] := sVal * V0 / dx;
            end;
          end;
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
    ApplySponge;   // после PEC, до записи пробника

    // ---------- ДОБАВЛЕНО: диагностика порта ----------
    // Ток в зазоре = ток смещения через ячейку зазора (стандартный приём FDTD).
    var gapE := Ex[idxEx(src_i, src_j - 1, src_k - 1)];
    var curI := eps0 * dy * dz * (gapE - prevGap) / dt;
    prevGap := gapE;

    WriteLn(gr, (n * dt):0:16, #9, gapE * dx:0:10, #9, curI:0:10);

    // Накопление ДПФ V и I на несущей — только в установившихся частях посылок
    // (после фронта и до начала спада), чтобы переходные процессы не портили Z.
    var inOn := false;
    for var p := 1 to nP do
      if (n >= onStart[p] + KeyRamp) and (n <= onEnd[p] - KeyRamp) then
        inOn := true;
    if inOn then
    begin
      ZreV += gapE * dx * Cos(w0 * n);
      ZimV -= gapE * dx * Sin(w0 * n);
      ZreI += curI  * Cos(w0 * n);
      ZimI -= curI  * Sin(w0 * n);
    end;

    // ---------- Запись сигнала пробника на каждом шаге ----------
    // ИЗМЕНЕНО: пробник сбоку от диполя, перпендикулярно оси — там излучение
    // максимально и поле ближе к волновому. (Было: src_i+7 вдоль оси.)
    WriteLn(pr, (n * dt):0:16, #9,
            Ex[idxEx(src_i, src_j + 5, src_k - 1)]:0:10);

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

  Close(pr);
  Close(gr);

  // ---------- ДОБАВЛЕНО: входной импеданс на несущей ----------
  // Z = V/I в комплексной арифметике: Z = V * conj(I) / |I|^2
  var den := Sqr(ZreI) + Sqr(ZimI);
  if den > 0 then
  begin
    var Zre := (ZreV * ZreI + ZimV * ZimI) / den;
    var Zim := (ZimV * ZreI - ZreV * ZimI) / den;

    WriteLn('Входной импеданс на f0: Z = ', Zre:0:1, ' + j*', Zim:0:1, ' Ом');
    WriteLn('  (мнимая ~ 0 -> резонанс; j<0 -> диполь коротковат, плечи +1;',
            ' j>0 -> длинноват)');
  end
  else
    WriteLn('Импеданс не посчитан: ток за посылки не накопился.');

  WriteLn('Расчёт завершён.');
end.
