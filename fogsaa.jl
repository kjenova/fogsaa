immutable AlignmentScoring
	match::Int
	mismatch::Int
	gap_start::Int
	gap::Int
end

immutable FOGSAABranch
	m::Int
	x::Int
	y::Int
end

function FOGSAA(X, Y, s::AlignmentScoring)
	const XY = 0x1
	const xY = 0x2
	const Xy = 0x3

	g0 = endof(X) - endof(Y)
	gM = g0 == 0? 0 : s.gap_start

	if g0 < 0
		gM += -g0 * s.gap
		m0 = endof(X)
	else
		gM += g0 * s.gap
		m0 = endof(Y)
	end

	gm = gM + m0 * s.mismatch
	gM += m0 * s.match

	alignments = zeros(UInt8, endof(X) + 1, endof(Y) + 1)
	scores = fill(gm, endof(X) + 1, endof(Y) + 1)
	scores[1, 1] = 0

	Mi = mi = gM - gm + 1
	queue = Vector{Vector{FOGSAABranch}}(Mi)
	queue[Mi] = [FOGSAABranch(1, 1, 1)]

	function enqueue(M, m, x, y, alignment)
		if alignment & 1 != 0
			x = -x
		end

		if alignment & 2 != 0
			y = -y
		end

		branch = FOGSAABranch(m, x, y)

		i = M - gm + 1

		if i > Mi
			Mi = i
		end

		if i < mi || mi == 0
			mi = i
		end

		if isdefined(queue, i)
			j = 1

			while j <= endof(queue[i]) && queue[i][j].m <= m
				j += 1
			end

			insert!(queue[i], j, branch)
		else
			queue[i] = [branch]
		end
	end

	function dequeue()
		x = queue[Mi][end].x
		y = queue[Mi][end].y

		pop!(queue[Mi])

		alignment = 0x0

		if x < 0
			x = -x
			alignment |= 1
		end

		if y < 0
			y = -y
			alignment |= 2
		end

		while !isdefined(queue, Mi) || endof(queue[Mi]) == 0
			if Mi == mi
				Mi = mi = 0
				break
			end

			Mi -= 1
		end

		x, y, alignment
	end

	ox = endof(X) + 1
	oy = endof(Y) + 1

	function print_alignment()
		X! = X[ox:end]
		Y! = Y[oy:end]			

		x = ox
		y = oy

		while x != 1 && y != 1
			if alignments[x, y] == xY
				y -= 1

				X! = string("-", X!)
				Y! = string(Y[y], Y!)
			elseif alignments[x, y] == Xy
				x -= 1

				X! = string(X[x], X!)
				Y! = string("-", Y!)
			else
				x -= 1
				y -= 1

				X! = string(X[x], X!)
				Y! = string(Y[y], Y!)
			end
		end

		if x != 1
			t = x - 1
			X! = string(X[1:t], X!)
			Y! = string(repeat("-", t), Y!)
		elseif y != 1
			t = y - 1
			Y! = string(Y[1:t], Y!)
			X! = string(repeat("-", t), X!)
		end

		r = endof(X!) - endof(Y!)

		if r < 0
			X! = string(X!, repeat("-", -r))
		else
			Y! = string(Y!, repeat("-", r))
		end

		println(X!)
		println(Y!)
	end

	while Mi - 1 + gm > scores[ox, oy]
	let
		x, y, alignment = dequeue()

		while x <= endof(X) && y <= endof(Y)
			sorted = 0x0
			score = [(scores[x, y] +
				  (X[x] == Y[y]? s.match : s.mismatch)),
				 (scores[x, y] + s.gap +
				  (alignment == xY? 0 : s.gap_start)),
				 (scores[x, y] + s.gap +
				  (alignment == Xy? 0 : s.gap_start))]
			M = Vector{Int}(3)
			m = Vector{Int}(3)

			if score[XY] > scores[x + 1, y + 1]
				scores[x + 1, y + 1] = M[XY] = score[XY]
			let
				rx = endof(X) - x
				ry = endof(Y) - y
				g = rx - ry
				r = ry

				if g != 0
					M[XY] += s.gap_start

					if g < 0
						g = -g
						r = rx
					end
				end

				M[XY] += g * s.gap
				m[XY] = M[XY] + r * s.mismatch
				M[XY] += r * s.match

				if M[XY] > scores[ox, oy]
					sorted = XY
				end
			end
			end

			if score[xY] > scores[x, y + 1]
				scores[x, y + 1] = M[xY] = score[xY]
			let
				rx = endof(X) + 1 - x
				ry = endof(Y) - y
				g = ry - rx
				r = rx

				if g < 0
					g = -g
					r = ry
					M[xY] += s.gap_start
				end

				M[xY] += g * s.gap
				m[xY] = M[xY] + r * s.mismatch
				M[xY] += r * s.match

				if M[xY] > scores[ox, oy]
					if sorted == 0x0 ||
					  (M[xY] >= M[XY] &&
					   (M[xY] > M[XY] ||
					    m[xY] > m[XY]))
						sorted <<= 2
						sorted |= xY
					else
						sorted |= xY << 2
					end
				end
			end
			end

			if score[Xy] > scores[x + 1, y]
				scores[x + 1, y] = M[Xy] = score[Xy]
			let
				rx = endof(X) - x
				ry = endof(Y) + 1 - y
				g = rx - ry
				r = ry

				if g < 0
					g = -g
					r = rx
					M[Xy] += s.gap_start
				end

				M[Xy] += g * s.gap
				m[Xy] = M[Xy] + r * s.mismatch
				M[Xy] += r * s.match

				if M[Xy] > scores[ox, oy]
					t = sorted >> 2

					if t == 0 ||
					  (M[Xy] >= M[t] &&
					   (M[Xy] > M[t] ||
					    m[Xy] > m[t]))
						u = sorted & 0b11

						if u == 0 ||
						  (M[Xy] >= M[u] &&
						   (M[Xy] > M[u] ||
						    m[Xy] > m[u]))
							sorted <<= 2
							sorted |= Xy
						else
							sorted = u
							t = t << 4 | Xy << 2
							sorted |= t
						end
					else
						sorted |= Xy << 4
					end
				end
			end
			end

			alignment = sorted & 0b11

			if alignment == 0x0
				@goto prune
			end

			a = (sorted & 0b1100) >> 2

			if a != 0 && M[a] > m[alignment]
				enqueue(M[a], m[a], a == xY? x : x + 1,
					a == Xy? y : y + 1, a)

				b = sorted >> 4

				if b != 0 && M[b] > m[a]
					enqueue(M[b], m[b], b == xY? x : x + 1,
						b == Xy? y : y + 1, b)
				end
			end

			x = alignment == xY? x : x + 1
			y = alignment == Xy? y : y + 1
			alignments[x, y] = alignment
		end

		ox = x
		oy = y
@label prune
	end
	end

	print_alignment()
end
