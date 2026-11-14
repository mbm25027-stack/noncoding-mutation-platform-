import pathlib
import re
import shutil
import sys

# Optional: paste your Google Cloud OAuth 2.0 Web Client ID here
GOOGLE_CLIENT_ID = "YOUR_GOOGLE_CLIENT_ID.apps.googleusercontent.com"
OWNER_EMAIL = "shiladityabiswas22@gmail.com"

target = pathlib.Path("index.html")
if not target.exists():
    alt = list(pathlib.Path("files (1)").glob("*.html")) if pathlib.Path("files (1)").exists() else []
    if alt:
        shutil.copy2(alt[0], target)
        print(f"Copied {alt[0]} -> index.html")
    else:
        sys.exit("Error: index.html not found in current folder.")

shutil.copy2(target, "index.html.bak")
html = target.read_text(encoding="utf-8")

if "fb-google-btn-slot" in html:
    print("index.html is already upgraded!")
    sys.exit(0)

# 1. Add Google Identity Services script in <head>
assert "</head>" in html, "Could not find </head>"
html = html.replace(
    "</head>",
    '<script src="https://accounts.google.com/gsi/client" async defer></script>\n</head>',
    1,
)

# 2. Add Google Sign-In styles before </style>
google_css = """.fb-google-row { margin-top: 16px; display: flex; flex-direction: column; gap: 10px; }
.fb-or { display: flex; align-items: center; gap: 10px; font: 500 10.5px var(--mono); letter-spacing: .12em; text-transform: uppercase; color: #7f93a8; }
.fb-or::before, .fb-or::after { content: ""; flex: 1; height: 1px; background: rgba(255,255,255,.14); }
.fb-g-btn { display: inline-flex; align-items: center; justify-content: center; gap: 10px; width: 100%; padding: 11px 16px; border-radius: 8px; border: 1px solid rgba(255,255,255,.22); background: #ffffff; color: #0f172a; font: 600 14px var(--sans); cursor: pointer; transition: transform .2s, box-shadow .2s; }
.fb-g-btn:hover { transform: translateY(-1px); box-shadow: 0 8px 20px rgba(0,0,0,.25); }
.fb-g-pill { display: none; align-items: center; justify-content: space-between; gap: 10px; padding: 10px 14px; border-radius: 8px; background: rgba(101,210,187,.12); border: 1px solid rgba(101,210,187,.4); color: #e6ecf8; font: 500 13px var(--sans); }
.fb-g-pill.on { display: flex; }
.fb-g-info { display: flex; align-items: center; gap: 10px; min-width: 0; }
.fb-g-info img { width: 28px; height: 28px; border-radius: 50%; border: 1px solid #65d2bb; }
.fb-g-info span { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.fb-g-out { background: none; border: 0; color: #65d2bb; font: 600 11.5px var(--mono); cursor: pointer; text-decoration: underline; flex-shrink: 0; }
</style>"""
assert "</style>" in html, "Could not find </style>"
html = html.replace("</style>", google_css, 1)

# 3. Upgrade Question 1 (data-step="0") with Google Sign-In
step0_pattern = re.compile(
    r'<div class="fb-slide active" data-step="0">.*?</div>\s*(?=<div class="fb-slide" data-step="1">)',
    re.DOTALL,
)
new_step0 = """<div class="fb-slide active" data-step="0">
            <div class="fb-q">Question 1</div>
            <h3><label for="fb-name">What&rsquo;s your name?</label></h3>
            <input class="fb-input" id="fb-name" type="text" autocomplete="name" maxlength="80" placeholder="Your full name">
            <div class="fb-google-row">
              <div class="fb-or">or sign in with Google</div>
              <div id="fb-google-btn-slot">
                <button type="button" class="fb-g-btn" id="fb-google-fallback">
                  <svg width="18" height="18" viewBox="0 0 48 48" aria-hidden="true"><path fill="#EA4335" d="M24 9.5c3.54 0 6.71 1.22 9.21 3.6l6.85-6.85C35.9 2.38 30.47 0 24 0 14.62 0 6.51 5.38 2.56 13.22l7.98 6.19C12.43 13.72 17.74 9.5 24 9.5z"/><path fill="#4285F4" d="M46.98 24.55c0-1.57-.15-3.09-.38-4.55H24v9.02h12.94c-.58 2.96-2.26 5.48-4.78 7.18l7.73 6c4.51-4.18 7.09-10.36 7.09-17.65z"/><path fill="#FBBC05" d="M10.53 28.59c-.48-1.45-.76-2.99-.76-4.59s.27-3.14.76-4.59l-7.98-6.19C.92 16.46 0 20.12 0 24c0 3.88.92 7.54 2.56 10.78l7.97-6.19z"/><path fill="#34A853" d="M24 48c6.48 0 11.93-2.13 15.89-5.81l-7.73-6c-2.15 1.45-4.92 2.3-8.16 2.3-6.26 0-11.57-4.22-13.47-9.91l-7.98 6.19C6.51 42.62 14.62 48 24 48z"/></svg>
                  Continue with Google
                </button>
              </div>
              <div class="fb-g-pill" id="fb-google-user">
                <div class="fb-g-info">
                  <img id="fb-g-avatar" src="" alt="">
                  <span id="fb-g-label"></span>
                </div>
                <button type="button" class="fb-g-out" id="fb-g-out">Sign out</button>
              </div>
            </div>
          </div>

          """
html, n0 = step0_pattern.subn(new_step0, html, count=1)
assert n0 == 1, "Could not match data-step='0' block"

# 4. Upgrade Thank You slide (data-step="4")
step4_pattern = re.compile(
    r'<div class="fb-slide fb-done" data-step="4">.*?</div>\s*(?=</div>\s*<div class="fb-err")',
    re.DOTALL,
)
new_step4 = f"""<div class="fb-slide fb-done" data-step="4">
            <svg class="fb-check" viewBox="0 0 84 84" aria-hidden="true"><circle cx="42" cy="42" r="30"/><path d="M29 43l9 9 18-20"/></svg>
            <h3>Thank you<span id="fb-thanks-name"></span>!</h3>
            <p id="fb-status-msg">Your feedback has been sent to <strong>{OWNER_EMAIL}</strong>.</p>
            <div class="fb-links"><a id="fb-mail" href="#">Open in email client</a><button type="button" id="fb-again">Send another</button></div>
          </div>
"""
html, n4 = step4_pattern.subn(new_step4, html, count=1)
assert n4 == 1, "Could not match data-step='4' block"

# 5. Upgrade Feedback JavaScript
fb_script_pattern = re.compile(r"<script>\s*/\* Feedback card-slider.*?</script>", re.DOTALL)
new_fb_script = f"""<script>
/* Feedback card-slider + Gmail delivery + Google Sign-In */
(function () {{
  var OWNER_EMAIL = '{OWNER_EMAIL}';
  var FEEDBACK_ENDPOINT = 'https://formsubmit.co/ajax/' + OWNER_EMAIL;
  var GOOGLE_CLIENT_ID = '{GOOGLE_CLIENT_ID}';

  var $ = function (id) {{ return document.getElementById(id); }};
  var track = $('fb-track'); if (!track) return;
  var slides = track.querySelectorAll('.fb-slide'), step = 0, TOTAL = 4;
  var googleProfile = null;

  var words = ['', 'Very poor', 'Poor', 'Below average', 'Fair', 'Okay', 'Good', 'Very good', 'Great', 'Excellent', 'Outstanding'];
  var rate = $('fb-rate'), h = '';
  for (var i = 1; i <= 10; i++) h += '<label class="fb-num"><input type="radio" name="fb-sat" value="' + i + '" aria-label="' + i + ' out of 10"><span>' + i + '</span></label>';
  rate.innerHTML = h;
  rate.addEventListener('change', function (e) {{ $('fb-word').textContent = words[+e.target.value] || ''; $('fb-err').textContent = ''; }});
  document.getElementById('fb-purpose-h').parentNode.addEventListener('change', function () {{ $('fb-err').textContent = ''; }});

  function decodeJwt(token) {{
    var base64Url = token.split('.')[1];
    var base64 = base64Url.replace(/-/g, '+').replace(/_/g, '/');
    var jsonPayload = decodeURIComponent(atob(base64).split('').map(function (c) {{
      return '%' + ('00' + c.charCodeAt(0).toString(16)).slice(-2);
    }}).join(''));
    return JSON.parse(jsonPayload);
  }}

  function applyGoogleUser(profile) {{
    googleProfile = profile;
    if (profile.name && !$('fb-name').value.trim()) $('fb-name').value = profile.name;
    $('fb-g-avatar').src = profile.picture || '';
    $('fb-g-label').textContent = profile.name + ' (' + profile.email + ')';
    $('fb-google-user').classList.add('on');
    $('fb-google-btn-slot').style.display = 'none';
    $('fb-err').textContent = '';
  }}

  function clearGoogleUser() {{
    googleProfile = null;
    $('fb-google-user').classList.remove('on');
    $('fb-google-btn-slot').style.display = '';
    if (window.google && google.accounts && google.accounts.id) google.accounts.id.disableAutoSelect();
  }}

  $('fb-g-out').addEventListener('click', clearGoogleUser);

  function initGoogleSignIn() {{
    var hasClientId = GOOGLE_CLIENT_ID && GOOGLE_CLIENT_ID.indexOf('YOUR_GOOGLE_CLIENT_ID') === -1;
    if (hasClientId && window.google && google.accounts && google.accounts.id) {{
      google.accounts.id.initialize({{
        client_id: GOOGLE_CLIENT_ID,
        callback: function (res) {{
          try {{
            var p = decodeJwt(res.credential);
            applyGoogleUser({{ name: p.name || '', email: p.email || '', picture: p.picture || '' }});
          }} catch (e) {{ fail('Could not verify Google account. Please enter your name above.'); }}
        }}
      }});
      var slot = $('fb-google-btn-slot');
      slot.innerHTML = '';
      google.accounts.id.renderButton(slot, {{ theme: 'outline', size: 'large', width: 320, text: 'continue_with' }});
    }} else {{
      var fbBtn = $('fb-google-fallback');
      if (fbBtn) {{
        fbBtn.addEventListener('click', function () {{
          if (window.google && google.accounts && google.accounts.id && hasClientId) {{
            google.accounts.id.prompt();
          }} else {{
            fail('Add your Google OAuth Client ID in index.html to activate Google Sign-In, or type your name above.');
          }}
        }});
      }}
    }}
  }}
  window.addEventListener('load', initGoogleSignIn);

  var reveal = [$('fb-intro'), $('fb-stack')];
  if ('IntersectionObserver' in window) {{
    var io = new IntersectionObserver(function (es) {{ es.forEach(function (e) {{ if (e.isIntersecting) {{ e.target.classList.add('in-view'); io.unobserve(e.target); }} }}); }}, {{ threshold: 0.2 }});
    reveal.forEach(function (el) {{ io.observe(el); }});
  }} else reveal.forEach(function (el) {{ el.classList.add('in-view'); }});

  function val(id) {{ return $(id).value.trim(); }}
  function checked(n) {{ var c = document.querySelector('input[name="' + n + '"]:checked'); return c ? c.value : ''; }}
  function setInert() {{ slides.forEach(function (sl, i) {{ var on = i === step; sl.classList.toggle('active', on); if (on) {{ sl.removeAttribute('inert'); sl.removeAttribute('aria-hidden'); }} else {{ sl.setAttribute('inert', ''); sl.setAttribute('aria-hidden', 'true'); }} }}); }}
  function go(n) {{
    step = n; setInert();
    track.style.transform = 'translateX(-' + (step * 100) + '%)';
    var shown = Math.min(step + 1, TOTAL);
    $('fb-count').textContent = '0' + shown + ' / 0' + TOTAL;
    $('fb-prog').style.width = (step >= TOTAL ? 100 : shown / TOTAL * 100) + '%';
    $('fb-back').disabled = step === 0 || step >= TOTAL;
    $('fb-nav').style.display = step >= TOTAL ? 'none' : '';
    $('fb-next').innerHTML = step === TOTAL - 1 ? 'Submit &#10003;' : 'Next &rarr;';
    $('fb-err').textContent = '';
    setTimeout(function () {{ var f = slides[step].querySelector('input'); if (f && step < TOTAL) f.focus({{ preventScroll: true }}); }}, 450);
  }}
  function fail(msg) {{ var e = $('fb-err'); e.textContent = msg; e.classList.remove('shake'); void e.offsetWidth; e.classList.add('shake'); }}
  function valid() {{
    if (step === 0 && !val('fb-name')) {{ fail('Please tell us your name or sign in with Google.'); return false; }}
    if (step === 1 && !checked('fb-purpose')) {{ fail('Please choose a purpose.'); return false; }}
    if (step === 3 && !checked('fb-sat')) {{ fail('Please pick a score from 1 to 10.'); return false; }}
    return true;
  }}
  function submit() {{
    var d = {{
      _subject: 'NC-MutScan Feedback from ' + val('fb-name'),
      _template: 'table',
      name: val('fb-name'),
      google_email: googleProfile ? googleProfile.email : '(not signed in with Google)',
      purpose: checked('fb-purpose'),
      organisation: val('fb-org') || '(not given)',
      satisfaction: checked('fb-sat') + '/10',
      submitted: new Date().toISOString()
    }};
    try {{ var all = JSON.parse(localStorage.getItem('ncmutscan_feedback') || '[]'); all.push(d); localStorage.setItem('ncmutscan_feedback', JSON.stringify(all)); }} catch (e) {{}}

    $('fb-status-msg').innerHTML = 'Sending feedback to <strong>' + OWNER_EMAIL + '</strong>&hellip;';
    fetch(FEEDBACK_ENDPOINT, {{
      method: 'POST',
      headers: {{ 'Content-Type': 'application/json', 'Accept': 'application/json' }},
      body: JSON.stringify(d)
    }}).then(function (r) {{
      if (r.ok) {{
        $('fb-status-msg').innerHTML = 'Delivered directly to <strong>' + OWNER_EMAIL + '</strong>. Thank you!';
      }} else {{
        $('fb-status-msg').innerHTML = 'Saved! Click <strong>Open in email client</strong> below to send directly to <strong>' + OWNER_EMAIL + '</strong>.';
      }}
    }}).catch(function () {{
      $('fb-status-msg').innerHTML = 'Saved! Click <strong>Open in email client</strong> below to send directly to <strong>' + OWNER_EMAIL + '</strong>.';
    }});

    var body = 'Name: ' + d.name + '\\nGoogle Email: ' + d.google_email + '\\nPurpose: ' + d.purpose + '\\nOrganisation: ' + d.organisation + '\\nSatisfaction: ' + d.satisfaction;
    $('fb-mail').href = 'mailto:' + OWNER_EMAIL + '?subject=' + encodeURIComponent(d._subject) + '&body=' + encodeURIComponent(body);
    $('fb-thanks-name').textContent = d.name ? ', ' + d.name.split(' ')[0] : '';
    go(TOTAL);
  }}
  $('fb-next').addEventListener('click', function () {{ if (!valid()) return; if (step === TOTAL - 1) submit(); else go(step + 1); }});
  $('fb-back').addEventListener('click', function () {{ if (step > 0) go(step - 1); }});
  $('fb-again').addEventListener('click', function () {{
    ['fb-name', 'fb-org'].forEach(function (id) {{ $(id).value = ''; }});
    document.querySelectorAll('#feedback input[type=radio]').forEach(function (r) {{ r.checked = false; }});
    $('fb-word').textContent = ''; go(0);
  }});
  ['fb-name', 'fb-org'].forEach(function (id) {{ $(id).addEventListener('keydown', function (e) {{ if (e.key === 'Enter') {{ e.preventDefault(); $('fb-next').click(); }} }}); }});
  setInert(); go(0);
}})();
</script>"""

html, ns = fb_script_pattern.subn(new_fb_script, html, count=1)
assert ns == 1, "Could not match Feedback <script> block"

target.write_text(html, encoding="utf-8")
print("SUCCESS: index.html updated (backup saved at index.html.bak).")