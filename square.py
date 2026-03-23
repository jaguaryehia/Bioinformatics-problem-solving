def create_square(row, col, piece=None):
    alphacols = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h'}
    return {'row': row, 'col': col, 'piece': piece, 'alphacol': alphacols[col]}

def eq_square(square1, square2):
    return square1['row'] == square2['row'] and square1['col'] == square2['col']

def has_piece(square):
    return square['piece'] is not None

def isempty(square):
    return not has_piece(square)

def has_team_piece(square, color):
    return has_piece(square) and square['piece'].color == color

def has_enemy_piece(square, color):
    return has_piece(square) and square['piece'].color != color

def isempty_or_enemy(square, color):
    return isempty(square) or has_enemy_piece(square, color)

def in_range(*args):
    for arg in args:
        if arg < 0 or arg > 7:
            return False

    return True

def get_alphacol(col):
    alphacols = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h'}
    return alphacols[col]
