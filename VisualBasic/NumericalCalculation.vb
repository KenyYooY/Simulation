Option Explicit

'0   1   3   -5  -9
'3   -2  1   1   6
'8   -5  -7  2   -15
'1   1   -2  -3  -15

Public Function GaussJordan(ByVal rngA As Range, Optional ByVal eps As Double = 0.000001) As Variant

    Dim A As Variant
    Dim aa As Double, b As Double, s As Double
    Dim n As Integer, n1 As Integer
    ' Dim X() As Double
    Dim X As Variant
    Dim ma As Integer
    Dim fmax As Double, fmax1 As Double
    Dim i As Integer, ii As Integer, i1 As Integer, j As Integer, l As Integer
    
    
    A = rngA
    If UBound(A, 2) <= UBound(A, 1) Then
        Error (0)
    End If
    
    n = UBound(A, 1)
    n1 = n + 1
    
    ' 行方向、列方向、両方の範囲選択での回答に対応するため2次元配列を用意
    ReDim X(n - 1, n - 1)
    
    For i = 1 To n
        ' Pivot
        ma = i
        fmax = Math.Abs(A(i, i))
        If i + 1 < n Then
            i1 = i + 1
            For ii = i1 To n
                fmax1 = Math.Abs(A(ii, i))
                If fmax1 > fmax Then
                    ma = ii
                    fmax = fmax1
                End If
            Next ii
        End If
        
        ' Error : Too small Diagonal
        If fmax <= eps Then
            Error (0)
        End If
        
        ' Exchenge
        If ma <> i Then
            For j = i To n1
                aa = A(ma, j)
                A(ma, j) = A(i, j)
                A(i, j) = aa
            Next j
        End If
        
        ' Solve
        b = A(i, i)
        
        For j = 1 To n1
            A(i, j) = A(i, j) / b
        Next j
        
        For l = 1 To n
            If l <> i Then
                s = A(l, i)
                For j = 1 To n1
                    A(l, j) = A(l, j) - s * A(i, j)
                Next j
            End If
        Next l
        
    Next i
    
    For i = 1 To n
        ' X(i - 1) = A(i, n1)
        ' 行方向、列方向、両方の範囲選択での回答に対応するため2次元配列に値を格納
        X(i - 1, 0) = A(i, n1)
        X(0, i - 1) = A(i, n1)
    Next i
    
    GaussJordan = X
    
End Function

