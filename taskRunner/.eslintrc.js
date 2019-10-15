module.exports = {
  parser: 'babel-eslint',
  extends: [
    'airbnb',
    'plugin:flowtype/recommended',
    'plugin:css-modules/recommended',
    'prettier',
    'prettier/flowtype',
    'prettier/react',
  ],
  plugins: ['flowtype', 'css-modules'],
  globals: { __DEV__: true },
  env: { browser: true },

  rules: {
    'max-len': [2, 120],
    'no-unused-vars': [1, { ignoreRestSiblings: true, argsIgnorePattern: '^_' }],
    'class-methods-use-this': 1,
    'import/first': 1,
    'no-plusplus': [2, { allowForLoopAfterthoughts: true }],
    'no-console': [1, { allow: ['warn', 'error', 'info'] }],
    'no-underscore-dangle': 0,
    'no-nested-ternary': 0,
    'import/named': 2,
    'react/prop-types': 0,
    // 'immutable/no-let': 2,
    // 'immutable/no-this': 2,
    // 'immutable/no-mutation': 2,
    'jsx-a11y/href-no-hash': 'off',
    'jsx-a11y/label-has-for': 1,

    'react/jsx-filename-extension': [2, { extensions: ['.js', '.jsx'] }],

    // Conversion is done automatically via
    'react/prefer-stateless-function': 'off',

    // 'prettier/prettier': [
    //   'warn',
    //   {
    //     singleQuote: true,
    //     trailingComma: 'all',
    //     printWidth: 120,
    //   },
    // ],
  },

  globals: {
    __DEV__: false,
  },
};
